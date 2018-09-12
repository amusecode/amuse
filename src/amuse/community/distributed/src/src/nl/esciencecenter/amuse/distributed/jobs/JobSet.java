/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nl.esciencecenter.amuse.distributed.jobs;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCreationFailedException;
import ibis.ipl.IbisFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Properties;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.pilots.PilotManager;
import nl.esciencecenter.amuse.distributed.pilots.PilotSet;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 */
public class JobSet extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(JobSet.class);

    public static final String PORT_NAME = "job.manager";

    private final Ibis ibis;

    //all pending jobs.
    private final LinkedList<AmuseJob> queue;

    //all pending/running/completed jobs.
    private final List<AmuseJob> jobs;

    private final PilotSet pilots;

    public JobSet(String serverAddress, PilotSet pilots) throws DistributedAmuseException {
        jobs = new ArrayList<AmuseJob>();

        this.pilots = pilots;

        logger.info("*** new jobset {}", serverAddress);

        try {
            Properties properties = new Properties();
            properties.put("ibis.server.address", serverAddress);
            properties.put("ibis.pool.name", "amuse");
            properties.put("ibis.location", "daemon@local");
            //properties.put("ibis.managementclient", "true");
            //properties.put("ibis.bytescount", "true");

            ibis = IbisFactory.createIbis(DistributedAmuse.IPL_CAPABILITIES, properties, true, pilots.getStatusMonitor(), null,
                    "master", DistributedAmuse.ONE_TO_ONE_PORT_TYPE, DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

            //label this ibis as the master node by running an election with us as the only 
            ibis.registry().elect("amuse");

            ibis.registry().enableEvents();

        } catch (IOException | IbisCreationFailedException e) {
            throw new DistributedAmuseException("failed to create ibis", e);
        }

        queue = new LinkedList<AmuseJob>();

        //start a thread to run the scheduling
        setName("Job Manager");
        setDaemon(true);
        start();
    }

    public Ibis getIbis() {
        return ibis;
    }

    public synchronized AmuseJob getJob(int jobID) throws DistributedAmuseException {
        for (AmuseJob job : jobs) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }

        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    //FUNCTION JOBS

    private synchronized void addJob(AmuseJob job) {
        queue.add(job);

        jobs.add(job);

        //run scheduler thread now
        notifyAll();
    }

    public FunctionJob submitFunctionJob(FunctionJobDescription description) throws DistributedAmuseException {
        FunctionJob result = new FunctionJob(description, ibis, this);

        addJob(result);

        return result;
    }

    //SCRIPT JOBS

    public ScriptJob submitScriptJob(ScriptJobDescription description) throws DistributedAmuseException {
        ScriptJob result = new ScriptJob(description, ibis, this);

        addJob(result);

        return result;

    }
    
    public WorkerJob submitWorkerJob(WorkerJobDescription jobDescription) throws DistributedAmuseException {
        WorkerJob result = new WorkerJob(jobDescription, ibis, this);

        addJob(result);

        return result;
    }

    public synchronized void removeJob(int jobID) throws DistributedAmuseException {
        //remove from queue
        for (int i = 0; i < queue.size(); i++) {
            if (queue.get(i).getJobID() == jobID) {
                queue.remove(i);
            }
        }
        
        for (int i = 0; i < jobs.size(); i++) {
            if (jobs.get(i).getJobID() == jobID) {
                AmuseJob result = jobs.remove(i);
                
                result.cancel();
                return;
            }
        }
        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    private synchronized boolean allScriptJobsDone() {
        for (AmuseJob job : jobs) {
            if (job.getType() == "script" && !job.isDone()) {
                return false;
            }
        }
        return true;
    }

    public synchronized void waitForScriptJobs() throws DistributedAmuseException {
        while (!allScriptJobsDone()) {
            try {
                wait();
            } catch (InterruptedException e) {
                throw new DistributedAmuseException("Interrupted while waiting for all jobs to finish");
            }
        }
    }

    public void end() {
        this.interrupt();

        for (AmuseJob job : getJobs()) {
            try {
                job.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel job: " + job, e);
            }
        }
        
        try {
            logger.debug("Terminating ibis pool");
            ibis.registry().terminate();
        } catch (IOException e) {
            logger.error("Failed to terminate ibis pool", e);
        }

        try {
            ibis.end();
        } catch (IOException e) {
            logger.error("Failed to end ibis", e);
        }
    }

    public synchronized AmuseJob[] getJobs() {
        return jobs.toArray(new AmuseJob[jobs.size()]);
    }

    /**
     * Wake up the scheduler thread.
     */
    public synchronized void nudge() {
        notifyAll();
    }

    /**
     * Scheduler thread. Periodically checks if suitable nodes can be found for jobs.
     */
    @Override
    public synchronized void run() {
        while (true) {
            int nodesInQueue = 0;
            //find nodes for jobs to run on
            Iterator<AmuseJob> iterator = queue.iterator();
            while (iterator.hasNext()) {
                AmuseJob job = iterator.next();

                if (job.isPending()) {
                    //find nodes to run this job on. Always only a single pilot, but may contain multiple nodes per pilot.
                    PilotManager target = pilots.getSuitablePilot(job);

                    //If suitable nodes are found
                    if (target != null) {
                        job.start(target);
                        //remove this job from the queue
                        iterator.remove();
                    } else {
                        nodesInQueue++;
                    }
                } else {
                    //remove this job from the queue
                    iterator.remove();
                }
            }

            if (nodesInQueue > 0) {
                logger.info("Now " + nodesInQueue + " waiting in queue");
            }
            
            try {
                wait(5000);
            } catch (InterruptedException e) {
                logger.debug("Scheduler thread interrupted, time to quit");
                return;
            }
            
            
        }
    }

    public synchronized List<Integer> getWorkerIDs() {
        ArrayList<Integer> result = new ArrayList<Integer>();
        
        for (AmuseJob job: getJobs()) {
            if (job.getType() == "worker") {
                result.add(job.getJobID());
            }
        }
        
        return result;
    }

    public FunctionJob getFunctionJob(int jobID) throws DistributedAmuseException {
        AmuseJob result = getJob(jobID);
        
        if (result.getType() != "function") {
            throw new DistributedAmuseException("Found job not a function job, instead found type " + result.getType()); 
        }
        
        return (FunctionJob) result;
    }

    public ScriptJob getScriptJob(int jobID) throws DistributedAmuseException {
        AmuseJob result = getJob(jobID);
        
        if (result.getType() != "script") {
            throw new DistributedAmuseException("Found job not a script job, instead found type " + result.getType()); 
        }
        
        return (ScriptJob) result;
    }
    
    public WorkerJob getWorkerJob(int jobID) throws DistributedAmuseException {
        AmuseJob result = getJob(jobID);
        
        if (result.getType() != "worker") {
            throw new DistributedAmuseException("Found job not a worker job, instead found type " + result.getType()); 
        }
        
        return (WorkerJob) result;
    }


}
