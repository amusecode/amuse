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

import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.WriteMessage;

import java.io.IOException;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;

/**
 * @author Niels Drost
 * 
 */
public class WorkerJob extends Job {

    private final WorkerDescription description;

    //only for worker jobs
    private ReceivePortIdentifier remoteWorkerPort = null;

    
    public WorkerJob(WorkerDescription description, Ibis ibis, JobManager jobManager) throws DistributedAmuseException {
        super(description.getNodeLabel(), description.getNrOfNodes(), ibis, jobManager);
        this.description = description;
    }

    public WorkerDescription getDescription() {
        return description;
    }
    
    private synchronized void setRemoteWorkerPort(ReceivePortIdentifier port) {
        this.remoteWorkerPort = port;
    }

    /**
     * @return
     */
    private synchronized IbisIdentifier[] getTargetIbisIdentifiers() {
        PilotNode [] target = getTarget();
        
        if (target == null) {
            return null;
        }

        IbisIdentifier[] result = new IbisIdentifier[target.length];

        for (int i = 0; i < result.length; i++) {
            result[i] = target[i].getIbisIdentifier();
        }

        return result;
    }

    @Override
    void writeJobDetails(WriteMessage writeMessage) throws IOException {
        writeMessage.writeObject(description);
        writeMessage.writeObject(getTargetIbisIdentifiers());
    }

    @Override
    void readJobStatus(ReadMessage readMessage) throws ClassNotFoundException, IOException {
//        ReceivePortIdentifier workerPort = (ReceivePortIdentifier) readMessage.readObject();
//        
//        setRemoteWorkerPort(workerPort);
    }

    /**
     * @param readMessage
     * @throws ClassNotFoundException
     * @throws IOException
     */
    @Override
    void readJobResult(ReadMessage readMessage) throws ClassNotFoundException, IOException {
        //NOTHING
    }



}
