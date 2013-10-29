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
package nl.esciencecenter.amuse.distributed.workers;

import ibis.ipl.Ibis;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.nio.ByteBuffer;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.JobManager;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class that handles incoming "worker" connections from AMUSE. Submits a job to the scheduler to allocate a node and start a
 * worker there, then forwards all messages to that node.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerConnectionServer extends Thread {

    public static final String WORKER_TYPE_STRING = "TYPE_WORKER";

    public static final String OUTPUT_TYPE_STRING = "TYPE_OUTPUT";

    public static final int TYPE_STRING_LENGTH = WORKER_TYPE_STRING.length();

    private static final Logger logger = LoggerFactory.getLogger(WorkerConnectionServer.class);

    private final ServerSocketChannel loopbackServer;

    private final WorkerOutputManager outputManager;

    private final Ibis ibis;

    private final JobManager scheduler;

    public WorkerConnectionServer(JobManager scheduler, File tmpDir) throws DistributedAmuseException {
        this.scheduler = scheduler;
        this.ibis = scheduler.getIbis();

        try {
            outputManager = new WorkerOutputManager(ibis);
            loopbackServer = ServerSocketChannel.open();
            // bind to a random port on local host
            loopbackServer.bind(new InetSocketAddress(InetAddress.getByName(null), 0));

            this.setName("worker connection server");
            this.setDaemon(true);
            this.start();
        } catch (IOException e) {
            throw new DistributedAmuseException("cannot start worker connection server", e);
        }
    }

    public int getPort() {
        return loopbackServer.socket().getLocalPort();
    }

    public void end() {
        try {
            loopbackServer.close();
        } catch (IOException e) {
            logger.error("Failed to close loopback server", e);
        }
    }

    public void run() {
        logger.debug("worker connection server started");
        while (true) {
            SocketChannel socket = null;
            try {
                socket = loopbackServer.accept();

                //turn on no-delay
                socket.socket().setTcpNoDelay(true);

                // read string, to make sure we are talking to amuse, and to get
                // the type of connection
                ByteBuffer magic = ByteBuffer.allocate(TYPE_STRING_LENGTH);

                while (magic.hasRemaining()) {
                    long read = socket.read(magic);

                    if (read == -1) {
                        throw new IOException("Connection closed on reading magic string");
                    }
                }

                String receivedString = new String(magic.array(), "UTF-8");
                if (receivedString.equalsIgnoreCase(WORKER_TYPE_STRING)) {
                    logger.debug("handling new worker connection");
                    new WorkerConnection(socket, ibis, scheduler);
                } else if (receivedString.equalsIgnoreCase(OUTPUT_TYPE_STRING)) {
                    logger.debug("handling new output connection");
                    outputManager.newOutputConnection(socket);
                } else {
                    throw new IOException("magic string (" + WORKER_TYPE_STRING + " or " + OUTPUT_TYPE_STRING
                            + ") not received. Instead got: " + receivedString);
                }

                logger.debug("New connection handled");
            } catch (Exception e) {
                if (socket != null) {
                    // report error to amuse
                    AmuseMessage errormessage = new AmuseMessage(0, AmuseMessage.FUNCTION_ID_INIT, 1, e.getMessage());
                    try {
                        errormessage.writeTo(socket);
                    } catch (Exception e2) {
                        // IGNORE
                    }
                    logger.error("error on starting remote code", e);

                    try {
                        socket.close();
                    } catch (Exception e2) {
                        // IGNORE
                    }
                }

                if (!loopbackServer.isOpen()) {
                    return;
                }

                // wait a bit before handling the next connection
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e1) {
                    // IGNORE
                }
            }
        }
    }

}
