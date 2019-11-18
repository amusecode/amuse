package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.ServerConnection;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.adaptors.schedulers.ssh.SshSchedulerAdaptor;
import nl.esciencecenter.xenon.utils.StreamForwarder;
import nl.esciencecenter.xenon.schedulers.JobDescription;
import nl.esciencecenter.xenon.schedulers.JobStatus;
import nl.esciencecenter.xenon.schedulers.Scheduler;
import nl.esciencecenter.xenon.schedulers.Streams;
import nl.esciencecenter.xenon.utils.JavaJobDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Hub extends Thread {

    private static final long TIMEOUT = 60000;

    private static final Logger logger = LoggerFactory.getLogger(Hub.class);

    private final String resourceName;

//    private final Job job;
    private final String job;

    private final ServerConnection serverConnection;

    private final String address;

    private final StreamForwarder forwarder;

    private final Scheduler scheduler;

    private static JavaJobDescription createJobDesciption(ResourceManager resource, String[] hubAddresses)
            throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        //~ result.setInteractive(true);

        AmuseConfiguration configuration = resource.getConfiguration();

        //classpath
        List<String> classpath = result.getJavaClasspath();
        classpath.add(configuration.getAmuseHome().getPath() + "/community/distributed/data/");
        classpath.add(configuration.getAmuseHome().getPath() + "/community/distributed/data/*");

        result.setJavaMain("ibis.ipl.server.Server");

        //arguments
        List<String> javaArguments = result.getJavaArguments();
        javaArguments.add("--remote");
        //javaArguments.add("--hub-only");
        javaArguments.add("--port");
        javaArguments.add("0");

        String hubs = null;
        for (String hub : hubAddresses) {
            if (hubs == null) {
                hubs = hub;
            } else {
                hubs = hubs + "," + hub;
            }
        }

        if (hubs != null) {
            javaArguments.add("--hub-addresses");
            javaArguments.add(hubs);
        }

        return result;
    }

    public Hub(ResourceManager resource, AmuseConfiguration config, String[] hubs ) throws DistributedAmuseException {

        try {
            resourceName = resource.getName();

            JobDescription jobDescription = createJobDesciption(resource, hubs);

            logger.debug("starting hub on {} with job description {} with arguments {}", resourceName, jobDescription,
                    jobDescription.getArguments());

            if (resource.isLocal()) {
                scheduler = Scheduler.create("local", null, null, null);

                jobDescription.setQueueName("unlimited");
                jobDescription.setMaxRuntime(0);
            } else if (resource.getHubQueueName() != null) {
                scheduler = resource.getScheduler();

                jobDescription.setQueueName(resource.getHubQueueName());
                jobDescription.setMaxRuntime(resource.getHubTimeMinutes());

            } else {

                scheduler = Scheduler.create("ssh", resource.getLocation(), resource.getCredential(), resource.getProperties());

                jobDescription.setQueueName("unlimited");
                jobDescription.setMaxRuntime(0);
            }
            logger.debug("starting hub on {} using scheduler {}", resourceName, scheduler);

            Streams streams = scheduler.submitInteractiveJob(jobDescription); // of batchjob?
 
            job = streams.getJobIdentifier();

            scheduler.waitUntilRunning(job, 0);
            
            logger.debug("started job " + job);

            forwarder = new StreamForwarder(streams.getStderr(), System.err);

            serverConnection = new ServerConnection(streams.getStdout(), streams.getStdin(), System.out, "Hub at "
                    + resource.getName() + ": ", TIMEOUT, null);

            address = serverConnection.getAddress();

            logger.debug("hub on {} has address {}", resource.getName(), address);

            setDaemon(true);
            start();
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start hub on " + resource.getName() + ": " + e, e);
        }
    }

    public void run() {
        try {
            JobStatus status = scheduler.getJobStatus(job);

            while (true) {
                JobStatus oldStatus = status;
                status = scheduler.getJobStatus(job);

                if (oldStatus.getState() != status.getState()) {
                    logger.info("Status of {} now {}", this, status);
                } else {
                    logger.debug("Status of {} now {}", this, status);
                }
                
                Thread.sleep(10000);
            }
        } catch (XenonException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InterruptedException e) {
            return;
        }

    }

    public String getAddress() {
        return address;
    }

    void stopHub() {
        serverConnection.closeConnection();
        
        //interrupt status thread
        interrupt();
    }

    @Override
    public String toString() {
        return "Hub [resourceName=" + resourceName + ", address=" + address + "]";
    }
}
