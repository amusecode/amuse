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
package nl.esciencecenter.amuse.distributed.remote;

import ibis.ipl.IbisIdentifier;
import ibis.ipl.RegistryEventHandler;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class Watchdog implements RegistryEventHandler {

    public static final long WATCHDOG_TIMEOUT = 60000;

    private static final Logger logger = LoggerFactory.getLogger(Watchdog.class);

    private long deadline;
    private boolean terminated;

    private IbisIdentifier master = null;

    Watchdog() {
        deadline = System.currentTimeMillis() + WATCHDOG_TIMEOUT;
        terminated = false;
    }

    @Override
    public synchronized void died(IbisIdentifier ibis) {
        logger.debug("{} died", ibis);
        left(ibis);
    }

    @Override
    public synchronized void electionResult(String election, IbisIdentifier winner) {
        logger.debug("{} won election {}", winner, election);
        if (election.equals("master") && winner != null) {
            master = winner;
        }
    }

    public synchronized void gotPing() {
        logger.debug("Got ping from master");

        if (logger.isDebugEnabled()) {
            long remaining = deadline - System.currentTimeMillis();
            logger.debug("Resetting deadline. There were {} ms remaining before the watchdog would have expired.", remaining);
        }

        //update deadline
        deadline = System.currentTimeMillis() + WATCHDOG_TIMEOUT;
    }

    @Override
    public synchronized void joined(IbisIdentifier arg0) {
        //IGNORE
    }

    @Override
    public synchronized void left(IbisIdentifier ibis) {
        logger.debug("{} left", ibis);
        if (master != null && ibis.equals(master)) {
            logger.debug("master left, terminating pilot");
            terminated = true;
            notifyAll();
        }
    }

    /**
     * 
     */
    @Override
    public synchronized void poolClosed() {
        //IGNORE
    }

    @Override
    public synchronized void poolTerminated(IbisIdentifier ibis) {
        logger.debug("pool terminated by {}", ibis);
        terminated = true;
        notifyAll();
    }
    
    @Override
    public void gotSignal(String arg0, IbisIdentifier arg1) {
        //IGNORE
    }


    /**
     * Waits until the pool is terminated, or the watchdog timer expires.
     */
    public synchronized boolean waitUntilExpired() {
        long now = System.currentTimeMillis();

        //reset deadline
        deadline = now + WATCHDOG_TIMEOUT;

        logger.debug("Watchdog started");

        //Continuously wait until the deadline has expired. Signal function above will update deadline.
        while (!terminated && now < deadline) {

            //number of ms until timer goes off.
            long timeout = deadline - now;

            try {
                logger.error("Waiting for {} ms", timeout);
                wait(timeout);
            } catch (InterruptedException e) {
                //IGNORE
            }

            now = System.currentTimeMillis();
        }

        if (!terminated) {
            logger.error("Watchdog expired!");
        } else {
            logger.info("Pool terminated");
        }

        return terminated;
    }


}
