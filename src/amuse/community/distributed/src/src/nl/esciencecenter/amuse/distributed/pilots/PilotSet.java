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

import java.io.File;
import java.util.ArrayList;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJob;
import nl.esciencecenter.amuse.distributed.resources.ResourceManager;
import nl.esciencecenter.amuse.distributed.resources.ResourceSet;
import nl.esciencecenter.xenon.Xenon;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Connection of pilots to run jobs. May still be in a queue, or already running.
 * 
 * @author Niels Drost
 * 
 */
public class PilotSet {

    private static final Logger logger = LoggerFactory.getLogger(PilotSet.class);

    private final XenonJobStatusMonitor jobStatusMonitor;

    private final ResourceSet resourceManager;

    private final PilotStatusMonitor statusMonitor;

    private final ArrayList<PilotManager> pilots;

    private final Xenon xenon;

    private final File tmpDir;

    private final boolean debug;

    public PilotSet(Xenon xenon, ResourceSet resourceManager, File tmpDir, boolean debug) throws DistributedAmuseException {
        this.xenon = xenon;
        this.resourceManager = resourceManager;
        this.tmpDir = tmpDir;
        this.debug = debug;
        pilots = new ArrayList<PilotManager>();
        jobStatusMonitor = new XenonJobStatusMonitor(this, xenon);
        this.statusMonitor = new PilotStatusMonitor(this);
    }

    public synchronized PilotManager newPilot(String resourceName, String queueName, int nodeCount, int timeMinutes, int slots,
            String nodeLabel, String options) throws DistributedAmuseException {
        logger.debug("reserving new nodes: resource name = " + resourceName + " queue name = " + queueName
                + " number of nodes = " + nodeCount + " time (in minutes) = " + timeMinutes + " node label = " + nodeLabel);

        ResourceManager resource = resourceManager.getResource(resourceName);

        PilotManager result = new PilotManager(resource, queueName, nodeCount, timeMinutes, slots, nodeLabel, options,
                resourceManager.getIplServerAddress(), resourceManager.getHubAddresses(), xenon, tmpDir, debug);

        pilots.add(result);
        
        //trigger an update in the job statuses
        jobStatusMonitor.nudge();

        return result;
    }

    public synchronized PilotManager getPilot(int reservationID) throws DistributedAmuseException {
        for (PilotManager reservation : pilots) {
            if (reservation.getAmuseID() == reservationID) {
                return reservation;
            }
        }
        throw new DistributedAmuseException("Reservation with ID " + reservationID + " not found");
    }

    public synchronized void deletePilot(int pilotID) throws DistributedAmuseException {
        logger.debug("deleting reservation " + pilotID);

        for (int i = 0; i < pilots.size(); i++) {
            PilotManager pilot = pilots.get(i);
            if (pilotID == pilot.getAmuseID()) {
                pilots.remove(i);
                pilot.stop();
                return;
            }
        }
        throw new DistributedAmuseException("Pilot " + pilotID + " not found");
    }

    private synchronized int[] getPilotIDs() {
        int[] result = new int[pilots.size()];

        for (int i = 0; i < result.length; i++) {
            result[i] = pilots.get(i).getAmuseID();
        }

        return result;
    }

    //check if all pilots are running. will throw an error if any pilots are done/failed
    private boolean allPilotsRunning() throws DistributedAmuseException {
        PilotManager[] pilots = getPilots();

        boolean result = true;

        for (PilotManager pilot : pilots) {
            if (pilot.hasException()) {
                throw new DistributedAmuseException("Pilot " + pilot + " failed while waiting for it to start",
                        pilot.getException());
            }
            if (pilot.isDone()) {
                throw new DistributedAmuseException("Pilot " + pilot + " done while waiting for it to start");
            }
            if (!pilot.isRunning()) {
                result = false;
            }
        }
        return result;
    }
    
    public synchronized void nudge() {
        notifyAll();
    }

    public void waitForAllPilots() throws DistributedAmuseException {
        logger.debug("waiting for all reservations to start");

        long timeout = 100; //ms

        int[] pilots = getPilotIDs();
        while (!allPilotsRunning()) {

            try {
                logger.debug("Now waiting {} ms", timeout);
                
                synchronized(this) {
                    wait(timeout);
                }
                //back-off waiting time to max 10 seconds
                if (timeout < 10000) {
                    timeout = timeout + 100;
                }
            } catch (InterruptedException e) {
                return;
            }
        }

        logger.debug("All {} pilots started.", pilots.length);
    }

    public synchronized void end() {
        for (PilotManager pilot : pilots) {
            try {
                pilot.stop();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to stop pilot: " + pilot, e);
            }
        }
    }

    public synchronized PilotManager[] getPilots() {
        return pilots.toArray(new PilotManager[pilots.size()]);
    }

    public synchronized PilotManager getPilot(UUID id) {
        for (PilotManager pilot : pilots) {
            if (pilot.getUniqueID().equals(id)) {
                return pilot;
            }
        }
        return null;
    }

    public PilotStatusMonitor getStatusMonitor() {
        return statusMonitor;
    }

    public PilotManager getSuitablePilot(AmuseJob job) {
        for (PilotManager pilot : getPilots()) {
            if (pilot.canRun(job)) {
                return pilot;
            }
        }
        return null;
    }
}
