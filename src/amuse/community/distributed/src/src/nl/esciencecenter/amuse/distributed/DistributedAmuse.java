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
package nl.esciencecenter.amuse.distributed;

import ibis.ipl.IbisCapabilities;
import ibis.ipl.PortType;

import java.io.File;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.jobs.JobSet;
import nl.esciencecenter.amuse.distributed.pilots.Lighthouse;
import nl.esciencecenter.amuse.distributed.pilots.PilotSet;
import nl.esciencecenter.amuse.distributed.resources.ResourceSet;
import nl.esciencecenter.amuse.distributed.web.WebInterface;
import nl.esciencecenter.amuse.distributed.workers.WorkerConnectionServer;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.amuse.distributed.util.Viz;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.qos.logback.classic.Level;

/**
 * Main Distributed AMUSE class. Started by AMUSE via the "Code" interface. Mostly contains objects that do the actual work.
 * 
 * @author Niels Drost
 * 
 */
public class DistributedAmuse {

    public static PortType MANY_TO_ONE_PORT_TYPE = new PortType(PortType.COMMUNICATION_RELIABLE, PortType.SERIALIZATION_OBJECT,
            PortType.RECEIVE_EXPLICIT, PortType.RECEIVE_AUTO_UPCALLS, PortType.CONNECTION_MANY_TO_ONE,
            PortType.CONNECTION_UPCALLS);

    public static PortType ONE_TO_ONE_PORT_TYPE = new PortType(PortType.COMMUNICATION_RELIABLE, PortType.SERIALIZATION_OBJECT,
            PortType.RECEIVE_EXPLICIT, PortType.RECEIVE_TIMEOUT, PortType.CONNECTION_ONE_TO_ONE);

    public static IbisCapabilities IPL_CAPABILITIES = new IbisCapabilities(IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED, IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS);

    private static final Logger logger = LoggerFactory.getLogger(DistributedAmuse.class);

    //resources potentially available for starting pilots on. Also starts hub on each resource, if required.
    private final ResourceSet resources;

    //pilots on resources available to run jobs 
    private final PilotSet pilots;

    //jobs running on some pilot running some worker or script.
    private final JobSet jobs;

    //talks to AMUSE, handling any requests to start workers.
    private final WorkerConnectionServer workerConnectionServer;

    //monitoring web interface
    private final WebInterface webInterface;

    private final boolean debug;

    private boolean ended = false;
    
    private final String amuseRootDir;
    
    private Viz viz;

    private static void initializeLogger(boolean debug) {
        if (debug) {
            ch.qos.logback.classic.Logger amuseLogger = (ch.qos.logback.classic.Logger) LoggerFactory
                    .getLogger("nl.esciencecenter.amuse");

            amuseLogger.setLevel(Level.DEBUG);

            logger.debug("DEBUG Enabled");
        }
    }

    public DistributedAmuse(UUID id, String codeDir, String amuseRootDir, int webInterfacePort, boolean debug, boolean startHubs,
            int workerQueueTimeout, int workerStartupTimeout) throws DistributedAmuseException {
        this.debug = debug;
        this.amuseRootDir = amuseRootDir;
        initializeLogger(debug);

        logger.info("Initializing Distributed Amuse {} with web interface on port {}, debug {}, and starting hubs {}",
                id.toString(), webInterfacePort, debug ? "enabled" : "disabled", startHubs ? "enabled" : "disabled");

        resources = new ResourceSet(amuseRootDir, startHubs);

        pilots = new PilotSet( resources,  id, debug);

        jobs = new JobSet(resources.getIplServerAddress(), pilots);
        
        workerConnectionServer = new WorkerConnectionServer(jobs, workerQueueTimeout, workerStartupTimeout);

        new Lighthouse(jobs.getIbis(), pilots);

        try {
            webInterface = new WebInterface(this, webInterfacePort);
        } catch (Exception e) {
            throw new DistributedAmuseException("could not create web interface", e);
        }
        logger.info("Distributed Amuse Initialized");
    }

    public boolean debugEnabled() {
        return debug;
    }

    public ResourceSet resources() {
        return resources;
    }

    public PilotSet pilots() {
        return pilots;
    }

    public JobSet jobs() {
        return jobs;
    }

    public WebInterface webInterface() {
        return webInterface;
    }

    /**
     * Port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     * 
     * @return the port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     */
    public int getWorkerPort() {
        logger.debug("Returning worker port.");
        return workerConnectionServer.getPort();
    }

    //to be on the safe side, synchronize this flag
    private synchronized boolean hasEnded() {
        return ended;
    }

    //to be on the safe side, synchronize this flag
    private synchronized void setEnded() {
        ended = true;
    }

    public synchronized void startup_viz() throws DistributedAmuseException {
        try {
          if ( viz != null) {
              viz.end();
              viz = null;
          }
          viz = new Viz(amuseRootDir, resources.getIplServerAddress());
        } catch (Exception e) {
            throw new DistributedAmuseException("could not create viz", e);
        }
    }

    /**
     * 
     */
    public void end() {
        if (hasEnded()) {
            return;
        }
        setEnded();
        
        try {
          if ( viz != null) {
              viz.end();
          }
        } catch (Exception e) {
            logger.info("could not shutdown viz", e);
        }

        logger.info("Ending distributed Amuse.");

        logger.debug("Ending web interface");
        webInterface.end();

        logger.debug("Ending worker connection server");
        workerConnectionServer.end();

        logger.debug("Ending job manager");
        jobs.end();

        //wait until all pilots have quit
        resources.endRegistry();

        logger.debug("Ending reservation manager");
        pilots.end();

        logger.debug("Ending resource manager");
        resources.end();

        logger.debug("Ending Xenon");

        logger.info("Distributed Amuse ended.");
    }

}
