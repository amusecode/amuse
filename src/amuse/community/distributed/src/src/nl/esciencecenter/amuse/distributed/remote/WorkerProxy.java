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

import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.lang.reflect.Field;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.SocketTimeoutException;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;
import java.util.Arrays;
import java.util.Map;

import javax.swing.event.ListSelectionEvent;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.WorkerJobDescription;
import nl.esciencecenter.amuse.distributed.workers.OutputForwarder;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * In charge of starting a AMUSE "worker" locally, and receiving requests from a WorkerConnection, and forwarding these to the
 * worker via a loopback socket.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerProxy extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(WorkerProxy.class);

    //how long until we give up on a worker initializing?
    private static final int ACCEPT_TIMEOUT = 100; // ms
    private static final int ACCEPT_TRIES = 50;

    static final String[] ENVIRONMENT_BLACKLIST = { "JOB_ID", "PE_", "PRUN_", "JOB_NAME", "JOB_SCRIPT", "OMPI_", "SLURM",
            "SBATCH", "SRUN", "HYDRA" };

    // local socket communication stuff

    private final SocketChannel socket;

    private final Process process;
    private int exitcode = 0;

    private final OutputForwarder out;
    private final OutputForwarder err;

    private final WorkerJobDescription description;
    private final AmuseConfiguration amuseConfiguration;

    //connection back to AMUSE

    private final Ibis ibis;

    private Exception error = null;

    private static Process startWorkerProcess(WorkerJobDescription description, AmuseConfiguration amuseConfiguration,
            int localSocketPort, File tempDirectory) throws Exception {
        File executable = new File(amuseConfiguration.getAmuseHome() + File.separator + description.getExecutable());

        if (!executable.canExecute()) {
            throw new DistributedAmuseException(executable + " is not executable, or does not exist");
        }

        ProcessBuilder builder = new ProcessBuilder();

        if (description.getNrOfThreads() > 0) {
            builder.environment().put("OMP_NUM_THREADS", Integer.toString(description.getNrOfThreads()));
        }

        if (!amuseConfiguration.isMPIEnabled()) {
            logger.info("not using mpiexec (as MPI is disabled in this AMUSE installation)");
            if (description.getNrOfWorkers() > 1) {
                throw new DistributedAmuseException("multiple workers (" + description.getNrOfWorkers()
                        + ") requested, but MPI disabled in this AMUSE installation");
            }

            // executable and port options
            builder.command().add(executable.toString());
            builder.command().add(Integer.toString(localSocketPort));

            //worker will connect to localhost
            builder.command().add("localhost");

            //mpi will _not_ be initialized in worker
            builder.command().add("false");
        } else {
            // use mpiexec
            builder.command().add(amuseConfiguration.getMpiexec());

            // set number of processes to start
            builder.command().add("-n");
            builder.command().add((Integer.toString(description.getNrOfWorkers())));

            // executable and port options
            builder.command().add(executable.toString());
            builder.command().add(Integer.toString(localSocketPort));

            //worker will connect to this machine (possibly from another machine)

            String localHostName = InetAddress.getLocalHost().getHostName();

            builder.command().add(localHostName);

            //mpi will be initialized in worker
            builder.command().add("true");
        }

        logger.info("starting worker process, command = " + builder.command());

        //start process and return
        return builder.start();
    }

    //small utility to figure out if the process is still running.
    private static boolean hasEnded(Process process) {
        try {
            process.exitValue();
            //we only end up here if the process is done
            return true;
        } catch (IllegalThreadStateException e) {
            //we got this exception as the process is not done yet
            return false;
        }
    }

    private static SocketChannel acceptConnection(ServerSocketChannel serverSocket, Process process) throws IOException,
            DistributedAmuseException {
        logger.debug("accepting connection");

        serverSocket.configureBlocking(true);
        serverSocket.socket().setSoTimeout(ACCEPT_TIMEOUT);

        for (int i = 0; i < ACCEPT_TRIES; i++) {
            try {

                //will timeout if this takes too long.
                //Note we do an accept on the socket, not the channel.
                //This is required as a workaround, otherwise the timeout on the accept is ignored
                SocketChannel result = serverSocket.socket().accept().getChannel();

                logger.debug("connection accepted");
                result.socket().setTcpNoDelay(true);
                return result;
            } catch (SocketTimeoutException e) {
                logger.debug("got timeout exception", e);
                if (hasEnded(process)) {
                    throw new DistributedAmuseException("worker failed to connect to java pilot, exited with exit code "
                            + process.exitValue());
                } else {
                    logger.debug("Got a timeout in accepting connection from worker. will keep trying");
                }
            }
        }
        throw new DistributedAmuseException("worker failed to connect to java pilot process within time");
    }

    /**
     * Starts a worker proxy. Make take a while.
     */
    public WorkerProxy(WorkerJobDescription description, AmuseConfiguration amuseConfiguration, Ibis ibis, File tempDirectory,
            int jobID) throws Exception {
        this.description = description;
        this.amuseConfiguration = amuseConfiguration;
        this.ibis = ibis;

        ServerSocketChannel serverSocket = ServerSocketChannel.open();
        //serverSocket.bind(new InetSocketAddress(InetAddress.getByName(null), 0), 10);
        serverSocket.bind(null);
        serverSocket.configureBlocking(true);
        //serverSocket.socket().setSoTimeout(ACCEPT_TIMEOUT);

        logger.debug("Bound server socket to " + serverSocket.socket().getLocalSocketAddress());

        //create process
        process = startWorkerProcess(description, amuseConfiguration, serverSocket.socket().getLocalPort(), tempDirectory);

        //attach streams
        out = new OutputForwarder(process.getInputStream(), description.getStdoutFile(), ibis);
        err = new OutputForwarder(process.getErrorStream(), description.getStderrFile(), ibis);

        logger.info("process started");

        socket = acceptConnection(serverSocket, process);
        serverSocket.close();

        logger.info("connection with local worker process established");

        //start a thread to start handling amuse requests
        setName("Worker Proxy for " + description.getID());
        setDaemon(true);
        start();
    }

    private synchronized void setError(Exception error) {
        this.error = error;
    }

    public synchronized Exception getError() {
        return error;
    }

    private synchronized void nativeKill() {
        if (process == null) {
            return;
        }

        try {
            Field f = process.getClass().getDeclaredField("pid");
            f.setAccessible(true);

            Object pid = f.get(process);

            ProcessBuilder builder = new ProcessBuilder("/bin/sh", "-c", "kill -9 " + pid.toString());
            
            builder.redirectError(Redirect.INHERIT);
            //builder.redirectInput();
            builder.redirectOutput(Redirect.INHERIT);

            logger.info("Killing process using command: " + Arrays.toString(builder.command().toArray()));
            
            Process killProcess = builder.start();
            
            killProcess.getOutputStream().close();
            
            int exitcode = killProcess.waitFor();
            
            logger.info("native kill done, result is " + exitcode);

        } catch (Throwable t) {
            logger.error("Error on (forcibly) killing process", t);
        }
    }

    public synchronized void end() {
        if (process != null) {
            process.destroy();

            try {
                exitcode = process.exitValue();
                logger.info("Process ended with result " + exitcode);
            } catch (IllegalThreadStateException e) {
                logger.error("Process not ended after process.destroy()! Trying native kill");
                nativeKill();
                try {
                    exitcode = process.exitValue();
                    logger.info("Process ended with result " + exitcode);
                } catch (IllegalThreadStateException e2) {
                    logger.error("Process not ended after native kill");
                }
            }
        }

        if (out != null) {
            // wait for out and err a bit
            try {
                out.join(1000);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }

        if (err != null) {
            try {
                err.join(1000);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }

    }

    /**
     * First connects to the "home" AMUSE ibis. Then continuously receives a message, performs a call, sends a reply.
     */
    public void run() {
        boolean running = true;
        long start, finish;

        AmuseMessage requestMessage = new AmuseMessage();
        AmuseMessage resultMessage = new AmuseMessage();

        SendPort sendPort = null;
        ReceivePort receivePort = null;

        try {
            sendPort = ibis.createSendPort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE);
            receivePort = ibis.createReceivePort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE, description.getID());

            receivePort.enableConnections();

            // get address of amuse node
            logger.debug("waiting for result of amuse election");
            IbisIdentifier amuse = ibis.registry().getElectionResult("amuse");

            //create a connection back to the amuse process via the ibis there.
            logger.debug("connecting to receive port of worker at amuse node");
            sendPort.connect(amuse, description.getID());
            logger.debug("connected, saying hello");
            WriteMessage helloMessage = sendPort.newMessage();
            helloMessage.writeObject(receivePort.identifier());
            helloMessage.writeObject(amuseConfiguration.getAmuseHome().getAbsolutePath().toString());
            helloMessage.finish();

            while (running) {
                logger.debug("Receiving call message");
                ReadMessage readMessage = receivePort.receive();

                if (readMessage == null) {
                    throw new IOException("cannot get request from worker");
                }

                logger.debug("Reading call request from IPL message");

                requestMessage.readFrom(readMessage);

                readMessage.finish();

                int functionID = requestMessage.getFunctionID();

                if (functionID == AmuseMessage.FUNCTION_ID_STOP) {
                    // final request handled
                    running = false;
                }

                // clean message
                resultMessage.clear();

                logger.debug("Performing call for function " + functionID);

                // perform call. Will put result in result message
                logger.debug("performing call with function ID " + requestMessage.getFunctionID());
                start = System.currentTimeMillis();
                requestMessage.writeTo(socket);
                resultMessage.readFrom(socket);
                finish = System.currentTimeMillis();

                logger.debug("done performing call with function ID " + requestMessage.getFunctionID() + " error = "
                        + resultMessage.getError());

                logger.debug("result: " + resultMessage);

                WriteMessage writeMessage = sendPort.newMessage();

                resultMessage.writeTo(writeMessage);

                writeMessage.writeLong(finish - start);

                writeMessage.finish();

                logger.debug("Done performing call for function " + functionID);

                if (logger.isDebugEnabled()) {
                    logger.debug("Call took " + (finish - start) + " ms");
                }

            }
        } catch (Exception e) {
            logger.error("Error while handling request, stopping worker", e);
            setError(e);
        }
        if (sendPort != null) {
            try {
                sendPort.close();
            } catch (IOException e) {
                //IGNORE
            }
        }
        if (receivePort != null) {
            try {
                receivePort.close();
            } catch (IOException e) {
                //IGNORE
            }
        }
        logger.debug("Worker proxy done");
    }

}
