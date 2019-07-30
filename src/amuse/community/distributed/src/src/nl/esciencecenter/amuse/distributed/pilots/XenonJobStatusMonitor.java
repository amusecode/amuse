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
package nl.esciencecenter.amuse.distributed.pilots;
import nl.esciencecenter.xenon.schedulers.JobStatus;
import nl.esciencecenter.xenon.XenonException;


/**
 * @author Niels Drost
 * 
 *         Utility class to track the status of one or more jobs. Will periodically poll the status of all pilot Xenon jobs, and
 *         set this status in the pilot. Requesting all statuses at once is much more efficient than one at a time.
 * 
 */
public class XenonJobStatusMonitor extends Thread {

    private final int JOB_UPDATE_DELAY = 5000;

    private final PilotSet pilotSet;

    XenonJobStatusMonitor(PilotSet pilotSet) {
        this.pilotSet = pilotSet;

        setName("Xenon Job Status monitor");
        setDaemon(true);
        start();
    }

    public synchronized void nudge() {
        notifyAll();
    }

    public void run() {
        while (true) {
            //Retrieve current set of pilots.
            PilotManager[] pilots = pilotSet.getPilots();

            //Create a list of xenon jobs.
            String[] jobs = new String[pilots.length];


            //Only add jobs for pilots which are still running.
            for (int i = 0; i < pilots.length; i++) {
                JobStatus oldStatus = pilots[i].getXenonJobStatus();

                if (oldStatus != null && oldStatus.isDone()) {
                    //this job is already done
                    jobs[i] = null;
                } else {
                    jobs[i] = pilots[i].getXenonJob();
                }
            }

            //Fetch all statuses from xenon
            JobStatus[] statuses = new JobStatus[ pilots.length];
            try {
                for (int i = 0; i < pilots.length; i++) {
                    if (jobs[i] != null) {
                        statuses[i]=pilots[i].getResource().getScheduler().getJobStatus(jobs[i]);
                    } else {
                        statuses[i]=null;
                    }
                }
            } catch (XenonException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            

            //Forward resulting statuses (if one is availabe) to the pilots
            for (int i = 0; i < pilots.length; i++) {
                if (statuses[i] != null) {
                    pilots[i].setXenonJobStatus(statuses[i]);
                }
            }
            
            //notify the pilot set, in case it is waiting for the pilots to start
            pilotSet.nudge();

            synchronized (this) {
                try {
                    wait(JOB_UPDATE_DELAY);
                } catch (InterruptedException e) {
                    //thread interrupted, end status monitor
                    return;
                }
            }
        }
    }

}
