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
package nl.esciencecenter.amuse.distributed.reservations;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.jobs.Job;
import nl.esciencecenter.xenon.jobs.JobStatus;

/**
 * @author Niels Drost
 * 
 *         Utility class to track the status of one or more jobs. Will periodically poll the status of jobs, and allow retrieval
 *         of the last-known status of jobs. Will also keep track of the last status of a job when it finishes.
 * 
 */
public class JobStatusMonitor extends Thread {

    private final int JOB_UPDATE_DELAY = 5000;

    private final Map<Job, JobStatus> statusMap;

    private final Xenon xenon;

    JobStatusMonitor(Xenon xenon) {
        this.xenon = xenon;
        statusMap = new HashMap<Job, JobStatus>();

        setName("Job Status monitor");
        setDaemon(true);
        start();
    }

    public synchronized void addJob(Job job) {
        statusMap.put(job, null);
        //trigger update to fill in status
        notifyAll();
    }

    public synchronized void removeJob(Job job) {
        statusMap.remove(job);
    }

    public synchronized JobStatus getstatus(Job job) {
        return statusMap.get(job);
    }

    private synchronized Job[] getInProgressJobs() {
        ArrayList<Job> result = new ArrayList<Job>();

        for (Map.Entry<Job, JobStatus> entry : statusMap.entrySet()) {
            if (entry.getValue() == null || !entry.getValue().isDone()) {
                result.add(entry.getKey());
            }
        }

        return result.toArray(new Job[result.size()]);
    }

    private synchronized void update(Job[] jobs, JobStatus[] statuses) {
        for (int i = 0; i < jobs.length; i++) {
            statusMap.put(jobs[i], statuses[i]);
        }
        notifyAll();
    }

    public synchronized void waitUntilAllRunning() throws DistributedAmuseException {
        while (true) {
            boolean allRunning = true;

            for (JobStatus status : statusMap.values()) {
                if (status == null || !status.isRunning()) {
                    allRunning = false;
                }
                if (status != null && status.hasException()) {
                    throw new DistributedAmuseException("Reservation failed while waiting for all reservations to start",
                            status.getException());
                }
                if (status != null && status.isDone()) {
                    throw new DistributedAmuseException("Reservation already done waiting for all reservations to start");
                }
            }

            if (allRunning) {
                return;
            }

            try {
                wait();
            } catch (InterruptedException e) {
                //return immediately if we get interrupted
                return;
            }
        }

    }

    public void run() {
        while (true) {
            Job[] jobs = getInProgressJobs();

            JobStatus[] statuses = xenon.jobs().getJobStatuses(jobs);

            update(jobs, statuses);

            try {
                Thread.sleep(JOB_UPDATE_DELAY);
            } catch (InterruptedException e) {
                //thread interrupted, end status monitor
                return;
            }
        }
    }

}
