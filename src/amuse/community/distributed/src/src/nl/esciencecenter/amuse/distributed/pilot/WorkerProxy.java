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
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.SocketTimeoutException;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.WorkerDescription;
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

    private static final String[] ENVIRONMENT_BLACKLIST = { "JOB_ID", "PE_", "PRUN_", "JOB_NAME", "JOB_SCRIPT", "OMPI_" };

    // local socket communication stuff

    private final SocketChannel socket;

    private final Process process;
    private int exitcode = 0;

    private final OutputForwarder out;
    private final OutputForwarder err;

    private final WorkerDescription description;
    private final AmuseConfiguration amuseConfiguration;

    //connection back to AMUSE

    private final Ibis ibis;

    private Exception error = null;

    /**
     * List of all hosts used, to give to MPI. Contains duplicates for all machines running multiple worker processes
     */
    private static String[] createHostnameList(WorkerDescription description, IbisIdentifier[] nodes)
            throws DistributedAmuseException {
        String[] hostnames = new String[description.getNrOfWorkers()];
        int next = 0;

        int nrOfProcesses = description.getNrOfWorkers();
        int nrOfNodes = description.getNrOfNodes();

        if (nrOfNodes != nodes.length) {
            throw new DistributedAmuseException("number of nodes used (" + nodes.length
                    + ") not equals to number of nodes required (" + nrOfNodes + ")");
        }

        hostnames = new String[nrOfProcesses];

        if (nodes.length == 1) {
            //just us, job done.
            for (int i = 0; i < nrOfProcesses; i++) {
                hostnames[i] = "localhost";
            }
        } else {
            for (int i = 0; i < nodes.length; i++) {
                String tags[] = nodes[i].tagAsString().split(",");
                if (tags.length != 3) {
                    throw new DistributedAmuseException("Cannot get tag from node identifier: " + nodes[i].tagAsString());
                }
                String hostname = tags[2];

                // number of processes per node
                int ppn = nrOfProcesses / nrOfNodes;
                // nrOfWorkers not dividable by nrOfNodes. see if this is
                // a "remainder" node with an extra worker
                if (i < nrOfProcesses % nrOfNodes) {
                    ppn++;
                }
                for (int j = 0; j < ppn; j++) {
                    hostnames[next] = hostname;
                    next++;
                }
            }
            if (next != nrOfProcesses) {
                logger.error("error in setting hostnames. List is of length " + next + " but should be " + nrOfProcesses);
            }
        }
        return hostnames;
    }

    private static File createHostFile(String[] hostnames, File tempDirectory) throws IOException {
        File result = File.createTempFile("host", "file", tempDirectory).getAbsoluteFile();

        FileWriter hostFileWriter = new FileWriter(result);
        for (String hostname : hostnames) {
            hostFileWriter.write(hostname + "\n");
        }
        hostFileWriter.flush();
        hostFileWriter.close();

        logger.info("host file = " + result);

        return result;
    }

    private static Process startWorkerProcess(WorkerDescription description, AmuseConfiguration amuseConfiguration,
            int localSocketPort, String[] hostnames, File tempDirectory) throws Exception {
        File executable = new File(amuseConfiguration.getAmuseHome() + File.separator + description.getExecutable());

        if (!executable.canExecute()) {
            throw new DistributedAmuseException(executable + " is not executable, or does not exist");
        }

        ProcessBuilder builder = new ProcessBuilder();

        File workingDirectory = new File(tempDirectory, "worker-" + description.getID());
        workingDirectory.mkdirs();
        builder.directory(workingDirectory);

        // make sure there is an "output" directory for a code to put output in
        //new File(workingDirectory, "output").mkdir();

        for (String key : builder.environment().keySet().toArray(new String[0])) {
            for (String blacklistedKey : ENVIRONMENT_BLACKLIST) {
                if (key.startsWith(blacklistedKey)) {
                    builder.environment().remove(key);
                    logger.info("removed " + key + " from environment");
                }
            }
        }
        if (description.getNrOfThreads() > 0) {
            builder.environment().put("OMP_NUM_THREADS", Integer.toString(description.getNrOfThreads()));
        }

        if (!amuseConfiguration.isMpiexecEnabled()) {
            logger.info("not using mpiexec (as it is disabled)");
            if (description.getNrOfWorkers() > 1) {
                throw new DistributedAmuseException("multiple workers (" + description.getNrOfWorkers()
                        + ") requested, but mpiexec disabled in this AMUSE installation");
            }
        } else if (description.getNrOfNodes() == 1) {
            // no need for machine file, set number of processes.
            builder.command(amuseConfiguration.getMpiexec(), "-n", Integer.toString(description.getNrOfWorkers()));
        } else {
            // use machine file
            File hostFile = createHostFile(hostnames, tempDirectory);
            builder.command(amuseConfiguration.getMpiexec(), "-machinefile", hostFile.getAbsolutePath());
        }

        // executable and port options
        builder.command().add(executable.toString());
        builder.command().add(Integer.toString(localSocketPort));

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
        } catch (IllegalStateException e) {
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
    public WorkerProxy(WorkerDescription description, AmuseConfiguration amuseConfiguration, IbisIdentifier[] nodes, Ibis ibis,
            File tempDirectory, int jobID) throws Exception {
        this.description = description;
        this.amuseConfiguration = amuseConfiguration;
        this.ibis = ibis;

        String[] hostnames = createHostnameList(description, nodes);

        ServerSocketChannel serverSocket = ServerSocketChannel.open();
        serverSocket.bind(new InetSocketAddress(InetAddress.getByName(null), 0));

        //create process
        process = startWorkerProcess(description, amuseConfiguration, serverSocket.socket().getLocalPort(), hostnames,
                tempDirectory);

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

    public synchronized void end() {
        if (process != null) {
            process.destroy();

            try {
                exitcode = process.exitValue();
                logger.info("Process ended with result " + exitcode);
            } catch (IllegalThreadStateException e) {
                logger.error("Process not ended after process.destroy()!");
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
                start = System.currentTimeMillis();
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
                requestMessage.writeTo(socket);
                resultMessage.readFrom(socket);
                logger.debug("done performing call with function ID " + requestMessage.getFunctionID() + " error = "
                        + resultMessage.getError());

                logger.debug("result: " + resultMessage);

                WriteMessage writeMessage = sendPort.newMessage();

                resultMessage.writeTo(writeMessage);

                writeMessage.finish();

                logger.debug("Done performing call for function " + functionID);
                finish = System.currentTimeMillis();

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
