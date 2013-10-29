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

import ibis.ipl.IbisIdentifier;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Node running one or more jobs. Note: NOT Thread-safe.
 * 
 * @author Niels Drost
 * 
 */
public class PilotNode {

    private static final Logger logger = LoggerFactory.getLogger(PilotNode.class);

    private final String label;

    private final int slots;

    private final int reservationID;

    private final String hostname;

    //address of this node
    private final IbisIdentifier ibisIdentifier;

    //list of all jobs running on this node
    private List<Job> jobs;

    public PilotNode(IbisIdentifier ibisIdentifier) throws DistributedAmuseException {
        this.ibisIdentifier = ibisIdentifier;

        jobs = new LinkedList<Job>();

        String[] tags = ibisIdentifier.tagAsString().split(",");

        if (tags.length != 4) {
            throw new DistributedAmuseException("Cannot parse ibis tag: " + ibisIdentifier.tagAsString());
        }

        try {
            reservationID = Integer.parseInt(tags[0]);
            label = tags[1];
            slots = Integer.parseInt(tags[2]);
            hostname = tags[3];

        } catch (NumberFormatException e) {
            throw new DistributedAmuseException("Cannot parse ibis tag: " + ibisIdentifier.tagAsString(), e);
        }
    }

    void addJob(Job job) {
        jobs.add(job);
    }

    boolean isAvailableForJobs() {
        int jobCount = 0;

        Iterator<Job> iterator = jobs.iterator();

        //ignore jobs that have already finished
        while (iterator.hasNext()) {
            Job job = iterator.next();

            if (job.isDone()) {
                iterator.remove();
            } else {
                jobCount++;
            }
        }

        logger.debug("pilot node {} currently running {} jobs, has {} slots available", this, jobCount, slots);
        
        return jobCount < slots;
    }
    
    public int getReservationID() {
        return reservationID;
    }

    public int getSlots() {
        return slots;
    }

    public String getLabel() {
        return label;
    }

    public IbisIdentifier getIbisIdentifier() {
        return ibisIdentifier;
    }

    public String getHostname() {
        return hostname;
    }

    @Override
    public String toString() {
        return "PilotNode [label=" + label + ", slots=" + slots + ", hostname=" + hostname + "]";
    }
}
