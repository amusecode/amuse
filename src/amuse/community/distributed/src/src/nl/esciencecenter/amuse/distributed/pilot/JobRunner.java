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
package nl.esciencecenter.amuse.distributed.pilot;

import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.jobs.WorkerDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A job running on a pilot node.S
 * 
 * @author Niels Drost
 * 
 */
public class JobRunner extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(JobRunner.class);

    private final int jobID;
    private final WorkerProxy workerProxy;
    private final ReceivePortIdentifier resultPort;
    private final Ibis ibis;

    /**
     * @param jobID
     * @param description
     * @param nodes
     * @param configuration
     * @param resultPort 
     * @throws Exception
     */
    public JobRunner(int jobID, WorkerDescription description, AmuseConfiguration configuration, IbisIdentifier[] nodes, ReceivePortIdentifier resultPort, Ibis ibis, File tmpDir)
            throws Exception {
        this.jobID = jobID;
        this.resultPort = resultPort;
        this.ibis = ibis;
        
        logger.debug("Starting job runner....");

        workerProxy = new WorkerProxy(description, configuration, nodes, ibis, tmpDir, jobID);

        setName("Job Runner for " + jobID);
    }

    public void run() {
        logger.debug("waiting for worker job {} to finish", jobID);

        try {
            workerProxy.join();
        } catch (InterruptedException e) {
            workerProxy.end();
        }

        logger.debug("worker {} done. Sending result to main amuse node.", jobID);

        //send result message to job
        try {
            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

            sendPort.connect(resultPort);

            WriteMessage message = sendPort.newMessage();

            message.writeObject(workerProxy.getError());
            message.finish();

            sendPort.close();

        } catch (IOException e) {
            logger.error("Failed to report status to main node", e);
        }

    }

}
