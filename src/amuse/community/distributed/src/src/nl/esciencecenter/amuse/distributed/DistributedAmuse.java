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

import nl.esciencecenter.amuse.distributed.jobs.JobManager;
import nl.esciencecenter.amuse.distributed.reservations.ReservationManager;
import nl.esciencecenter.amuse.distributed.resources.ResourceManager;
import nl.esciencecenter.amuse.distributed.web.WebInterface;
import nl.esciencecenter.amuse.distributed.workers.WorkerConnectionServer;
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.XenonFactory;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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

    //resources potentially available for starting reservations on. Also starts hub on each resource, if required.
    private final ResourceManager resourceManager;

    //starts pilots on resources. 
    private final ReservationManager reservationManager;

    //takes care of job queue, communicates with remote pilots
    private final JobManager jobManager;

    //talks to AMUSE, handling any worker requests and messages
    private final WorkerConnectionServer workerConnectionServer;

    //monitoring web interface
    private final WebInterface webInterface;

    //used to copy files, start jobs, etc.
    private final Xenon xenon;

    private File tmpDir;

    private static File createTmpDir() throws DistributedAmuseException {
        File systemTmpDir = new File(System.getProperty("java.io.tmpdir"));

        if (!systemTmpDir.exists()) {
            throw new DistributedAmuseException("Java tmpdir does not exist " + systemTmpDir);
        }

        File result = new File(systemTmpDir, "distributed-amuse/daemon-" + UUID.randomUUID().toString());
        result.mkdirs();

        return result;
    }

    public DistributedAmuse(String codeDir, String amuseRootDir, int webInterfacePort) throws DistributedAmuseException {
        logger.info("Initializing Distributed Amuse");
        try {
            xenon = XenonFactory.newXenon(null);
        } catch (XenonException  e) {
            throw new DistributedAmuseException("could not create Xenon library object", e);
        }

        tmpDir = createTmpDir();

        resourceManager = new ResourceManager(xenon, tmpDir, amuseRootDir);

        jobManager = new JobManager(resourceManager.getIplServerAddress(), tmpDir);
        
        reservationManager = new ReservationManager(xenon, resourceManager, jobManager.getNodes(), tmpDir);

        workerConnectionServer = new WorkerConnectionServer(jobManager, tmpDir);

        try {
            webInterface = new WebInterface(this, webInterfacePort);
        } catch (Exception e) {
            throw new DistributedAmuseException("could not create web interface", e);
        }
        logger.info("Distributed Amuse Initialized");
    }

    public ResourceManager resourceManager() {
        return resourceManager;
    }

    public ReservationManager reservationManager() {
        return reservationManager;
    }

    public JobManager jobManager() {
        return jobManager;
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

    /**
     * 
     */
    public void end() {
        logger.info("Ending distributed Amuse.");
        
        logger.debug("Ending web interface");
        webInterface.end();
        
        logger.debug("Ending worker connection server");
        workerConnectionServer.end();
        
        logger.debug("Ending job manager");
        jobManager.end();

        logger.debug("Ending registry");
        resourceManager.endRegistry();
        
        logger.debug("Ending reservation manager");
        reservationManager.end();
        
        logger.debug("Ending resource manager");
        resourceManager.end();
        
        logger.debug("Ending Xenon");
        XenonFactory.endAll();
        logger.info("Distributed Amuse ended.");
    }
    
}
