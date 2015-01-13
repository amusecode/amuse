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
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.SocketTimeoutException;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJobDescription;
import nl.esciencecenter.amuse.distributed.jobs.WorkerJobDescription;
import nl.esciencecenter.amuse.distributed.util.FileTransfers;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * In charge of starting a AMUSE "worker" locally, and receiving requests from a WorkerConnection, and forwarding these to the
 * worker via a loopback socket.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerJobRunner extends JobRunner {

    static final Logger logger = LoggerFactory.getLogger(WorkerJobRunner.class);

    //how long until we give up on a worker initializing?
    private static final int ACCEPT_TIMEOUT = 100; // ms

    static final String[] ENVIRONMENT_BLACKLIST = { "JOB_ID", "PE_", "PRUN_", "JOB_NAME", "JOB_SCRIPT", "OMPI_", "SLURM",
            "SBATCH", "SRUN", "HYDRA" };

    // local socket communication stuff

    private final SocketChannel socket;

    private final WorkerJobDescription description;

    //connection back to AMUSE

    private ProcessBuilder createProcessBuilder(int localSocketPort) throws Exception {
        Path executable;

        if (description.isDynamicPythonCode()) {
            //executable is the name of the worker without any directories.
            executable = sandbox.resolve(Paths.get(description.getExecutable()).getFileName());
            
            if (!Files.exists(executable)) {
                throw new DistributedAmuseException(executable + " does not exist");
            }
            
        } else {
            executable = Paths.get(description.getExecutable());

            //make absolute by prepending with amuse home
            if (!executable.isAbsolute()) {
                executable = amuseConfiguration.getAmuseHome().toPath().resolve(executable);
            }

            if (!Files.isExecutable(executable)) {
                throw new DistributedAmuseException(executable + " is not executable, or does not exist");
            }

        }
        
        ProcessBuilder builder = new ProcessBuilder();

        builder.directory(sandbox.toFile());

        if (description.getNrOfThreads() > 0) {
            builder.environment().put("OMP_NUM_THREADS", Integer.toString(description.getNrOfThreads()));
        }

        if (!amuseConfiguration.isMPIEnabled()) {
            logger.info("not using mpiexec (as MPI is disabled in this AMUSE installation)");
            if (description.getNrOfWorkers() > 1) {
                throw new DistributedAmuseException("multiple workers (" + description.getNrOfWorkers()
                        + ") requested, but MPI disabled in this AMUSE installation");
            }

            if (description.isDynamicPythonCode()) {
                Path amuseScriptPath = amuseConfiguration.getAmuseHome().getAbsoluteFile().toPath().resolve("amuse.sh");
                builder.command().add(amuseScriptPath.toString());
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

            if (description.isDynamicPythonCode()) {
                Path amuseScriptPath = amuseConfiguration.getAmuseHome().getAbsoluteFile().toPath().resolve("amuse.sh");
                builder.command().add(amuseScriptPath.toString());
            }
            
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
        return builder;
    }

    private SocketChannel acceptConnection(ServerSocketChannel serverSocket) throws IOException, DistributedAmuseException {
        logger.debug("accepting connection");

        serverSocket.configureBlocking(true);
        serverSocket.socket().setSoTimeout(ACCEPT_TIMEOUT);

        int accept_tries = (description.getStartupTimeout() * 1000) / ACCEPT_TIMEOUT;
        
        logger.debug("Will try to accept connection from worker " + accept_tries + " times");
        
        for (int i = 0; i < accept_tries; i++) {
            try {

                //will timeout if this takes too long.
                //Note we do an accept on the socket, not the channel.
                //This is required as a workaround, otherwise the timeout on the accept is ignored
                logger.debug("Waiting for connection");
                SocketChannel result = serverSocket.socket().accept().getChannel();

                logger.debug("connection accepted");
                result.socket().setTcpNoDelay(true);
                return result;
            } catch (SocketTimeoutException e) {
                logger.debug("got timeout exception", e);
                if (hasEnded()) {
                    throw new DistributedAmuseException("worker failed to connect to java pilot, exited with exit code "
                            + getExitCode());
                } else {
                    logger.debug("Got a timeout in accepting connection from worker. will keep trying");
                }
            }
        }
        throw new DistributedAmuseException("worker failed to connect to java pilot process within " + description.getStartupTimeout() + " seconds");
    }

    /**
     * Starts a worker proxy. Make take a while.
     */
    public WorkerJobRunner(AmuseJobDescription description, AmuseConfiguration configuration, ReceivePortIdentifier resultPort,
            Ibis ibis, Path tmpDir, boolean debug, ReadMessage message) throws Exception {
        super(description, configuration, resultPort, ibis, tmpDir, debug);

        this.description = (WorkerJobDescription) description;

        if (this.description.isDynamicPythonCode()) {
            FileTransfers.readDirectory(sandbox, message);
        }

        message.finish();

        
        InetSocketAddress wildcardAddress = new InetSocketAddress(0);
        
        ServerSocketChannel serverSocket = ServerSocketChannel.open();
        //serverSocket.bind(new InetSocketAddress(InetAddress.getByName(null), 0), 10);
        //serverSocket.bind(null);

        serverSocket.bind(wildcardAddress);
        serverSocket.configureBlocking(true);
        //serverSocket.socket().setSoTimeout(ACCEPT_TIMEOUT);

        logger.debug("Bound server socket to " + serverSocket.socket().getLocalSocketAddress());

        //create process
        startProcess(createProcessBuilder(serverSocket.socket().getLocalPort()));

        socket = acceptConnection(serverSocket);
        serverSocket.close();

        logger.info("connection with local worker process established");

        //start a thread to start handling amuse requests
        setName("Worker Job Runner for job " + description.getID());
        setDaemon(true);
        start();
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
            receivePort = ibis.createReceivePort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE, Integer.toString(description.getID()));

            receivePort.enableConnections();

            // get address of amuse node
            logger.debug("waiting for result of amuse election");
            IbisIdentifier amuse = ibis.registry().getElectionResult("amuse");

            //create a connection back to the amuse process via the ibis there.
            logger.debug("connecting to receive port of worker at amuse node");
            sendPort.connect(amuse, Integer.toString(description.getID()));
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
        logger.debug("Closing sendport");
        if (sendPort != null) {
            try {
                sendPort.close();
            } catch (IOException e) {
                //IGNORE
            }
        }
        logger.debug("Closing receiveport");
        if (receivePort != null) {
            try {
                receivePort.close();
            } catch (IOException e) {
                //IGNORE
            }
        }

        logger.debug("Waiting for process to end");
        waitForProcess();

        logger.debug("Worker job done, sending result to amuse");

        sendResult();
        
        logger.debug("Deleting sandbox");
        
        deleteSandbox();
    }

    @Override
    void writeResultData(WriteMessage message) throws IOException {
        //NOTHING
    }

}
