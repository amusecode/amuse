package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.ServerConnection;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.adaptors.ssh.SshAdaptor;
import nl.esciencecenter.xenon.credentials.Credential;
import nl.esciencecenter.xenon.engine.util.StreamForwarder;
import nl.esciencecenter.xenon.jobs.Job;
import nl.esciencecenter.xenon.jobs.JobDescription;
import nl.esciencecenter.xenon.jobs.JobStatus;
import nl.esciencecenter.xenon.jobs.Scheduler;
import nl.esciencecenter.xenon.jobs.Streams;
import nl.esciencecenter.xenon.util.JavaJobDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Hub extends Thread {

    private static final long TIMEOUT = 60000;

    private static final Logger logger = LoggerFactory.getLogger(Hub.class);

    private final String resourceName;

    private final Job job;

    private final ServerConnection serverConnection;

    private final String address;

    private final Xenon xenon;

    private static JavaJobDescription createJobDesciption(ResourceManager resource, String[] hubAddresses)
            throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        result.setInteractive(true);

        AmuseConfiguration configuration = resource.getConfiguration();

        //classpath
        List<String> classpath = result.getJavaClasspath();
        classpath.add(configuration.getAmuseHome().getPath() + "/src/amuse/community/distributed/src/dist/*");
        classpath.add(configuration.getAmuseHome().getPath() + "/src/amuse/community/distributed/worker.jar");
        classpath.add(configuration.getAmuseHome().getPath() + "/src/amuse/community/distributed");

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

    public Hub(ResourceManager resource, AmuseConfiguration config, String[] hubs, Xenon xenon) throws DistributedAmuseException {
        this.xenon = xenon;

        try {
            resourceName = resource.getName();

            JobDescription jobDescription = createJobDesciption(resource, hubs);

            logger.debug("starting hub on {} with job description {} with arguments {}", resourceName, jobDescription,
                    jobDescription.getArguments());

            Scheduler scheduler;

            if (resource.isLocal()) {
                scheduler = xenon.jobs().newScheduler("local", null, null, null);

                jobDescription.setQueueName("unlimited");
                jobDescription.setMaxTime(0);
            } else if (resource.getHubQueueName() != null) {
                scheduler = resource.getScheduler();

                jobDescription.setQueueName(resource.getHubQueueName());
                jobDescription.setMaxTime(resource.getHubTimeMinutes());

            } else {

                Credential credential = xenon.credentials().getDefaultCredential("ssh");

                Map<String, String> properties = new HashMap<String, String>();
                String gateway = resource.getGateway();
                if (gateway != null && !gateway.isEmpty()) {
                    properties.put(SshAdaptor.GATEWAY, gateway);
                }

                scheduler = xenon.jobs().newScheduler("ssh", resource.getLocation(), credential, properties);

                jobDescription.setQueueName("unlimited");
                jobDescription.setMaxTime(0);
            }
            logger.debug("starting hub on {} using scheduler {}", resourceName, scheduler);

            job = xenon.jobs().submitJob(scheduler, jobDescription);

            logger.debug("started job " + job);

            xenon.jobs().waitUntilRunning(job, 0);

            Streams streams = xenon.jobs().getStreams(job);

            new StreamForwarder(streams.getStderr(), System.err);

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
            JobStatus status = xenon.jobs().getJobStatus(job);

            while (true) {
                JobStatus oldStatus = status;
                status = xenon.jobs().getJobStatus(job);

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
