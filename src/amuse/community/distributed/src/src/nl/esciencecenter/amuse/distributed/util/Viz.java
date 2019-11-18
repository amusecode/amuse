package nl.esciencecenter.amuse.distributed.util;

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

public class Viz extends Thread {

    private static final long TIMEOUT = 60000;

    private static final Logger logger = LoggerFactory.getLogger(Viz.class);
    
    private final String address;

    private final String job;

    private final Scheduler scheduler; 

    public Viz(String amuseRootDir,String address ) throws DistributedAmuseException {

        try {
          
            this.address=address;
            JavaJobDescription jobDescription = new JavaJobDescription();

            //classpath
            List<String> classpath = jobDescription.getJavaClasspath();
            classpath.add(amuseRootDir + "/community/distributed/data/");
            classpath.add(amuseRootDir + "/community/distributed/data/*");

            jobDescription.setJavaMain("ibis.smartsockets.viz.SmartsocketsViz");

            //arguments
            List<String> javaArguments = jobDescription.getJavaArguments();
            javaArguments.add(address);

            logger.debug("starting Viz locally with job description {} with arguments {}", jobDescription,
                    jobDescription.getArguments());

            jobDescription.setQueueName("unlimited");
            jobDescription.setMaxRuntime(0);

            scheduler = Scheduler.create("local", null, null, null);

            logger.debug("starting Viz..");

            job = scheduler.submitBatchJob(jobDescription); // of batchjob?
 
            scheduler.waitUntilRunning(job, 0);
            
            logger.debug("started Viz as job " + job);

            setDaemon(true);
            start();
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start Viz " + ": " + e, e);
        }
    }

    public void end() throws DistributedAmuseException {
        logger.debug("cancelling xenon job for: {}", this);
        try {
            scheduler.cancelJob(job);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to cancel Viz job " + job, e);
        }
        this.stop();
    }

    @Override
    public String toString() {
        return "Viz [job=" + job + "]";
    }
}
