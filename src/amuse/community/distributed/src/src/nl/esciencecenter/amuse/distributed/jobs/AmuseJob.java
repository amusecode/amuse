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
import ibis.ipl.MessageUpcall;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.LinkedHashMap;
import java.util.Map;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.pilots.PilotManager;
import nl.esciencecenter.amuse.distributed.remote.Pilot;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A job run by the Distributed Amuse system. Contains description and status info, communicates with nodes that actually run job.
 * 
 * @author Niels Drost
 * 
 */
public abstract class AmuseJob extends Thread implements MessageUpcall {

    private enum State {
        PENDING, INITIALIZING, RUNNING, DONE, FAILED;
    }

    //how long we keep jobs which have finished.
    public static int TIMEOUT = 60000; //ms
    
    public static int RPC_TIMEOUT = 60000; //ms

    private static final Logger logger = LoggerFactory.getLogger(AmuseJob.class);

    private static int nextID = 0;

    static synchronized int getNextID() {
        return nextID++;
    }

    private final Ibis ibis;

    private final JobSet jobManager;

    private final ReceivePort resultReceivePort;

    private final AmuseJobDescription description;
    
    private State state;

    private PilotManager target = null;

    private Exception error = null;

    //will never timeout until timeout set
    private long expirationDate = Long.MAX_VALUE;

    public AmuseJob(AmuseJobDescription description, Ibis ibis, JobSet jobManager) throws DistributedAmuseException {
        this.description = description;

        this.ibis = ibis;
        this.jobManager = jobManager;

        this.state = State.PENDING;

        try {
            resultReceivePort = ibis
                    .createReceivePort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE, "job-" + description.getID(), this);
            resultReceivePort.enableConnections();
            resultReceivePort.enableMessageUpcalls();
        } catch (IOException e) {
            throw new DistributedAmuseException("Cannot create receive port for job", e);
        }
    }
    
    public String getType() {
        return description.getType();
    }

    public void end() {
        try {
            resultReceivePort.close();
        } catch (IOException e) {
            //IGNORE
        }
    }

    public AmuseJobDescription getDescription() {
        return description;
    }

    public int getJobID() {
        return description.getID();
    }

    private synchronized void setState(State newState) {
        state = newState;
        notifyAll();

        if (isDone()) {
            expirationDate = System.currentTimeMillis() + TIMEOUT;
        }
    }

    private synchronized void setError(Exception error) {
        this.error = error;
        setState(State.FAILED);
    }

    public synchronized String getJobState() {
        return state.toString();
    }

    public synchronized boolean isPending() {
        return state == State.PENDING;
    }

    public synchronized boolean isRunning() {
        return state == State.RUNNING;
    }

    public synchronized boolean isDone() {
        return state == State.DONE || state == State.FAILED;
    }

    public synchronized boolean hasFailed() {
        return state == State.FAILED;
    }

    public synchronized void waitUntilRunning(long timeout) throws Exception {
        long deadline;

        if (timeout == 0) {
            //wait for ever
            deadline = Long.MAX_VALUE;
        } else {
            deadline = System.currentTimeMillis() + timeout;
        }

        while (!isRunning()) {
            if (hasFailed()) {
                throw getError();
            }
            if (isDone()) {
                throw new DistributedAmuseException("Job already done while waiting until it is running");
            }
            long remaining = deadline - System.currentTimeMillis();
            if (remaining <= 0) {
                return;
            } else {
                try {
                    logger.debug("Waiting " + remaining + " millis for job to start running");
                    wait(remaining);
                } catch (InterruptedException e) {
                    return;
                }
            }
        }
    }

    public synchronized void waitUntilDone() {
        while (!isDone()) {
            logger.info("Waiting for {} job {} to finish, state now {}", getType(), getJobID(), getJobState());
            try {
                wait(5000);
            } catch (InterruptedException e) {
                return;
            }
        }
    }

    protected synchronized Exception getError() {
        return error;
    }

    protected synchronized PilotManager getTarget() {
        return target;
    }

    /**
     * @param target
     *            the nodes to run this job on.
     */
    public synchronized void start(PilotManager target) {
        if (!isPending()) {
            logger.error("Tried to run job {} that was not pending. Ignoring", this);
            return;
        }

        logger.debug("Running job {} on pilot {} with label {}", getJobID(), target.getID(), target.getLabel());

        //set state to initializing
        setState(State.INITIALIZING);

        this.target = target;

        target.addAmuseJob(this);

        //send out messages to the nodes in a separate thread (see run function below)
        setName("Job " + description.getID() + " starting thread");
        setDaemon(true);
        start();
    }

    /**
     * Function that starts the job.
     */
    @Override
    public void run() {
        try {
            PilotManager pilot = target;

            logger.debug("sending start command for {} to pilot", this);

            if (!pilot.isRunning()) {
                setError(new DistributedAmuseException("Pilot no longer running, cannot start job"));
                return;
            }

            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);
            ReceivePort receivePort = ibis.createReceivePort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE, null);
            receivePort.enableConnections();

            sendPort.connect(pilot.getIbisIdentifier(), "pilot", RPC_TIMEOUT, true);

            WriteMessage writeMessage = sendPort.newMessage();

            //command
            writeMessage.writeString("start");

            //where to send the reply for this message
            writeMessage.writeObject(receivePort.identifier());
            writeMessage.writeObject(this.resultReceivePort.identifier());
            writeMessage.writeObject(description);

            writeJobData(writeMessage);

            writeMessage.finish();

            sendPort.close();

            logger.debug("receiving reply from pilot");

            //FIXME: we should use some kind of rpc mechanism
            ReadMessage readMessage = receivePort.receive(RPC_TIMEOUT);

            Exception error = (Exception) readMessage.readObject();

            readMessage.finish();
            receivePort.close();

            if (error != null) {
                setError(new DistributedAmuseException("Remote node reported error:" + error, error));
                logger.debug("job {} ERROR", this);
            } else {
                setState(State.RUNNING);
                logger.debug("job {} started", this);
            }

        } catch (IOException | ClassNotFoundException e) {
            logger.error("Job failed!", e);
            setError(e);
        }
    }

    /**
     * Cancel a job by sending a message to the pilot
     */
    public void cancel() throws DistributedAmuseException {
        if (isDone()) {
            logger.trace("NOT Cancelling job as it is already done {}", this);
            return;
        } else if (isPending()) {
            logger.trace("Cancelling pending job {}", this);
            setError(new Exception("Job cancelled while pending"));
            return;
        } else {
            logger.debug("Cancelling job {}", this);
        }

        try {
            PilotManager master = target;

            if (!master.isRunning()) {
                setError(new DistributedAmuseException("Pilot no longer running, cannot cancel job"));
                return;
            }

            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

            sendPort.connect(master.getIbisIdentifier(), Pilot.PORT_NAME, RPC_TIMEOUT, true);

            WriteMessage writeMessage = sendPort.newMessage();

            //command
            writeMessage.writeString("cancel");

            writeMessage.writeInt(description.getID());

            writeMessage.finish();

            sendPort.close();

        } catch (IOException e) {
            throw new DistributedAmuseException("Failed to cancel job " + this, e);
        }
    }

    /**
     * Handles incoming result message from Pilots
     */
    @Override
    public void upcall(ReadMessage message) {
        logger.debug("Reading result message");

        Exception error;
        try {
            error = (Exception) message.readObject();

            //implemented by job sub type
            readJobResult(message);

        } catch (IOException | ClassNotFoundException e) {
            logger.error("Error while reading result", e);
            error = e;
        }

        if (error != null) {
            setError(new DistributedAmuseException("Remote node reported error: " + error, error));
        } else {
            setState(State.DONE);
        }

        jobManager.nudge();

        logger.debug("Result received, state now: {}", this);

    }

    public Map<String, String> getStatusMap() {
        Map<String, String> result = new LinkedHashMap<String, String>();

        result.put("ID", Integer.toString(description.getID()));
        result.put("Label", description.getLabel());
        result.put("State", getJobState().toString());
        result.put("Target", target.toString());

        Exception error = getError();

        if (error != null) {
            StringWriter writer = new StringWriter();
            PrintWriter printWriter = new PrintWriter(writer);
            error.printStackTrace(printWriter);
            result.put("Error", writer.toString());
        }

        return result;
    }

    abstract void writeJobData(WriteMessage writeMessage) throws IOException;

    abstract void readJobResult(ReadMessage readMessage) throws ClassNotFoundException, IOException;
    
    @Override
    public String toString() {
        return "AmuseJob [jobManager=" + jobManager + ", description=" + description + ", state=" + state + ", target=" + target
                + ", error=" + error + ", expirationDate=" + expirationDate + "]";
    }

}
