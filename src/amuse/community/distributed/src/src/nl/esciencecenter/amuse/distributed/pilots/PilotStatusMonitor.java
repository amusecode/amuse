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

import ibis.ipl.IbisIdentifier;
import ibis.ipl.RegistryEventHandler;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Handles events from the Ibis registry to update
 * 
 * @author Niels Drost
 * 
 */
public class PilotStatusMonitor implements RegistryEventHandler {

    private static final Logger logger = LoggerFactory.getLogger(PilotStatusMonitor.class);

    private final PilotSet pilots;

    public PilotStatusMonitor(PilotSet pilots) {
        this.pilots = pilots;
    }

    @Override
    public synchronized void died(IbisIdentifier ibis) {
        //handle like it left
        left(ibis);
    }

    @Override
    public void joined(IbisIdentifier ibis) {
        logger.debug("new Ibis joined: " + ibis);

        if (ibis.location().toString().equals("daemon@local")) {
            //ingore local deamon process
            return;
        }
        int id = Integer.parseInt(ibis.tagAsString());

        try {
            PilotManager pilot = pilots.getPilot(id);
            pilot.setIbisIdentifier(ibis);
        } catch (DistributedAmuseException e) {
            logger.error("Could not find matching PilotManager for joining Pilot " + id);
        }

        //wake up PilotSet to re-check if all pilots are now running
        pilots.nudge();
    }

    @Override
    public void left(IbisIdentifier ibis) {
        if (ibis.location().toString().equals("daemon@local")) {
            //ingore local deamon process
            return;
        }
        int id = Integer.parseInt(ibis.tagAsString());

        try {
            PilotManager pilot = pilots.getPilot(id);
            pilot.setLeft();
        } catch (DistributedAmuseException e) {
            logger.warn("Could not find matching PilotManager for leaving Pilot " + id);
        }
    }

    @Override
    public void electionResult(String name, IbisIdentifier winner) {
        //IGNORED
    }

    @Override
    public void gotSignal(String signal, IbisIdentifier origin) {
        //IGNORED
    }

    @Override
    public void poolClosed() {
        //IGNORED
    }

    @Override
    public void poolTerminated(IbisIdentifier arg0) {
        //IGNORED
    }
}
