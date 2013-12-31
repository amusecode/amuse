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

import java.io.File;
import java.util.ArrayList;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.PilotNodes;
import nl.esciencecenter.amuse.distributed.resources.Resource;
import nl.esciencecenter.amuse.distributed.resources.ResourceManager;
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.jobs.JobStatus;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Reservations of actual resources to run jobs. May still be in a queue, or already running.
 * 
 * @author Niels Drost
 * 
 */
public class ReservationManager {

    private final int JOIN_WAIT_TIME = 60000;

    private static final Logger logger = LoggerFactory.getLogger(ReservationManager.class);

    private final JobStatusMonitor jobStatusMonitor;

    private final ResourceManager resourceManager;

    private final PilotNodes nodes;

    private final ArrayList<Reservation> reservations;

    private final Xenon xenon;

    private final File tmpDir;

    public ReservationManager(Xenon xenon, ResourceManager resourceManager, PilotNodes nodes, File tmpDir)
            throws DistributedAmuseException {
        this.xenon = xenon;
        this.resourceManager = resourceManager;
        this.nodes = nodes;
        this.tmpDir = tmpDir;
        reservations = new ArrayList<Reservation>();
        jobStatusMonitor = new JobStatusMonitor(xenon);
    }

    public synchronized Reservation newReservation(String resourceName, String queueName, int nodeCount, int timeMinutes,
            int slots, String nodeLabel, String options) throws DistributedAmuseException {
        logger.debug("reserving new nodes: resource name = " + resourceName + " queue name = " + queueName
                + " number of nodes = " + nodeCount + " time (in minutes) = " + timeMinutes + " node label = " + nodeLabel);

        Resource resource = resourceManager.getResource(resourceName);

        Reservation result = new Reservation(resource, queueName, nodeCount, timeMinutes, slots, nodeLabel, options,
                resourceManager.getIplServerAddress(), resourceManager.getHubAddresses(), xenon, tmpDir);

        reservations.add(result);
        jobStatusMonitor.addJob(result.getJob());

        return result;
    }

    public synchronized Reservation getReservation(int reservationID) throws DistributedAmuseException {
        for (Reservation reservation : reservations) {
            if (reservation.getID() == reservationID) {
                return reservation;
            }
        }
        throw new DistributedAmuseException("Reservation with ID " + reservationID + " not found");
    }

    public synchronized void deleteReservation(int reservationID) throws DistributedAmuseException {
        logger.debug("deleting reservation " + reservationID);

        for (int i = 0; i < reservations.size(); i++) {
            Reservation reservation = reservations.get(i);
            if (reservationID == reservation.getID()) {
                reservations.remove(i);
                jobStatusMonitor.removeJob(reservation.getJob());
                reservation.cancel();
                return;
            }
        }
        throw new DistributedAmuseException("Reservation " + reservationID + " not found");
    }

    private synchronized int[] getReservationIDs() {
        int[] result = new int[reservations.size()];

        for (int i = 0; i < result.length; i++) {
            result[i] = reservations.get(i).getID();
        }

        return result;
    }

    public void waitForAllReservations() throws DistributedAmuseException {
        logger.debug("waiting for all reservations to start");

        long deadline = Long.MAX_VALUE;
        long timeout = 100; //ms
        
        int[] reservations = getReservationIDs();
        while (!nodes.containsNodesFrom(reservations)) {
            
            if (System.currentTimeMillis() > deadline) {
                throw new DistributedAmuseException("Pilot nodes failed to join");
            }
            //check for errors
            if (jobStatusMonitor.areAllRunning()) {
                //all jobs are now running. We only wait a fixed amount of time
                //until all the nodes are actually joined.
                if (deadline == Long.MAX_VALUE) {
                    logger.debug("All reservations running. Now waiting for all pilots to join");
                    deadline =  System.currentTimeMillis() + JOIN_WAIT_TIME;
                }
            }

            try {
                logger.debug("Now waiting {} ms", timeout);
                Thread.sleep(timeout);
                //back-off waiting time to max 10 seconds
                if (timeout < 10000) {
                    timeout = timeout * 2;
                }
            } catch (InterruptedException e) {
                //IGNORE
            }
        }

        logger.debug("All reservations started, all {} pilots accounted for.", reservations.length);
    }

    public synchronized void end() {
        for (Reservation reservation : reservations) {
            try {
                reservation.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel reservation: " + reservation, e);
            }
        }
    }

    public synchronized Reservation[] getReservations() {
        return reservations.toArray(new Reservation[reservations.size()]);
    }

    public String getState(Reservation reservation) {
        JobStatus status = jobStatusMonitor.getstatus(reservation.getJob());

        if (status == null) {
            return "UNKNOWN";
        } else {
            return status.getState();
        }
    }
}
