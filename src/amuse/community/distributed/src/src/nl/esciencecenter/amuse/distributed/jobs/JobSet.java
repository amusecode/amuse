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

    private final List<WorkerJob> workers;
    private final List<ScriptJob> scriptJobs;
    private final List<FunctionJob> functionJobs;

    private final PilotSet pilots;

    public JobSet(String serverAddress, PilotSet pilots, File tmpDir) throws DistributedAmuseException {
        workers = new ArrayList<WorkerJob>();
        scriptJobs = new ArrayList<ScriptJob>();
        functionJobs = new ArrayList<FunctionJob>();

        this.pilots = pilots;

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
        for (AmuseJob job : workers) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }

        for (AmuseJob job : functionJobs) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }

        for (AmuseJob job : scriptJobs) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }

        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    //FUNCTION JOBS

    private synchronized void addFunctionJob(FunctionJob job) {
        queue.add(job);

        functionJobs.add(job);

        //run scheduler thread now
        notifyAll();
    }

    public FunctionJob submitFunctionJob(FunctionJobDescription description) throws DistributedAmuseException {
        FunctionJob result = new FunctionJob(description, ibis, this);

        addFunctionJob(result);

        return result;
    }

    public synchronized FunctionJob[] getFunctionJobs() {
        return functionJobs.toArray(new FunctionJob[0]);
    }

    public synchronized FunctionJob getFunctionJob(int jobID) throws DistributedAmuseException {
        for (FunctionJob job : functionJobs) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }
        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    public void removeFunctionJob(int jobID) throws DistributedAmuseException {
        for (int i = 0; i < functionJobs.size(); i++) {
            if (functionJobs.get(i).getJobID() == jobID) {
                functionJobs.remove(i);
                return;
            }
        }
        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    //SCRIPT JOBS

    private synchronized void addScriptJob(ScriptJob job) {
        queue.add(job);

        scriptJobs.add(job);

        //run scheduler thread now
        notifyAll();
    }

    public ScriptJob submitScriptJob(ScriptJobDescription description) throws DistributedAmuseException {
        ScriptJob result = new ScriptJob(description, ibis, this);

        addScriptJob(result);

        return result;

    }

    public synchronized ScriptJob[] getScriptJobs() {
        return scriptJobs.toArray(new ScriptJob[0]);
    }

    public synchronized ScriptJob getScriptJob(int jobID) throws DistributedAmuseException {
        for (ScriptJob job : scriptJobs) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }
        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    public synchronized void removeScriptJob(int jobID) throws DistributedAmuseException {
        for (int i = 0; i < scriptJobs.size(); i++) {
            if (scriptJobs.get(i).getJobID() == jobID) {
                scriptJobs.remove(i);
                return;
            }
        }
        throw new DistributedAmuseException("Unknown job: " + jobID);
    }

    private synchronized boolean allScriptJobsDone() {
        for (ScriptJob job : scriptJobs) {
            if (!job.isDone()) {
                return false;
            }
        }
        return true;
    }

    //WORKER JOBS

    private synchronized void addWorkerJob(WorkerJob job) {
        queue.add(job);

        workers.add(job);

        //run scheduler thread now
        notifyAll();
    }

    public WorkerJob submitWorkerJob(WorkerJobDescription jobDescription) throws DistributedAmuseException {
        WorkerJob result = new WorkerJob(jobDescription, ibis, this);

        addWorkerJob(result);

        return result;
    }

    public synchronized WorkerJob[] getWorkerJobs() {
        return workers.toArray(new WorkerJob[0]);
    }

    public synchronized WorkerJob getWorkerJob(int jobID) throws DistributedAmuseException {
        for (WorkerJob job : workers) {
            if (jobID == job.getJobID()) {
                return job;
            }
        }
        throw new DistributedAmuseException("Unknown job: " + jobID);
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

        for (AmuseJob job : getWorkerJobs()) {
            try {
                job.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel job: " + job, e);
            }
        }
        for (AmuseJob job : getScriptJobs()) {
            try {
                job.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel job: " + job, e);
            }
        }
        for (AmuseJob job : getFunctionJobs()) {
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
                    }
                } else {
                    //remove this job from the queue
                    iterator.remove();
                }

            }

            try {
                wait(5000);
            } catch (InterruptedException e) {
                logger.debug("Scheduler thread interrupted, time to quit");
                return;
            }
        }
    }

    public synchronized int getWorkerJobCount() {
        return workers.size();
    }

    public synchronized int[] getWorkerIDs() {
        int[] result = new int[workers.size()];

        for (int i = 0; i < result.length; i++) {
            result[i] = workers.get(i).getJobID();
        }
        return result;
    }

}
