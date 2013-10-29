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
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Takes care of getting all output to its final destination. Either a local file, or the AMUSE process.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerOutputManager extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(WorkerOutputManager.class);

    public static final int MAX_MESSAGE_SIZE = 10240;

    private final Ibis ibis;
    private final ReceivePort receivePort;

    private final Map<UUID, SocketChannel> outputConnections;

    WorkerOutputManager(Ibis ibis) throws IOException {
        this.ibis = ibis;

        outputConnections = new HashMap<UUID, SocketChannel>();

        this.receivePort = ibis.createReceivePort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE, "output");
        this.receivePort.enableConnections();

        setDaemon(true);
        start();
    }

    private synchronized void removeOutputConnection(UUID id) {
        outputConnections.remove(id);
    }

    private synchronized SocketChannel getOutputConnection(UUID id) {
        return outputConnections.get(id);
    }

    void newOutputConnection(SocketChannel channel) throws IOException {
        UUID newID = UUID.randomUUID();

        synchronized (this) {
            outputConnections.put(newID, channel);
        }

        // send id of this connection to AMUSE
        AmuseMessage uuidReply = new AmuseMessage();

        uuidReply.addString(newID.toString());

        uuidReply.writeTo(channel);
    }

    @Override
    public void run() {
        ByteBuffer buffer = ByteBuffer.allocateDirect(MAX_MESSAGE_SIZE);

        while (!ibis.registry().hasTerminated()) {
            try {
                ReadMessage message = receivePort.receive();

                String file = message.readString();

                int count = message.readInt();

                buffer.clear();
                buffer.limit(count);

                message.readByteBuffer(buffer);

                message.finish();

                UUID fileID = null;

                try {
                    fileID = UUID.fromString(file);
                } catch (IllegalArgumentException e) {
                    // Apparently, this is not a UUID at all
                }

                buffer.flip();

                if (fileID == null) {
                    logger.trace("Got {} bytes for file {}", count, file);
                    // not a UUID, just write it to a file
                    FileOutputStream out = new FileOutputStream(file, true);
                    out.getChannel().write(buffer);
                    out.flush();
                    out.close();
                } else {
                    logger.trace("Got {} bytes for stream {}" , count, fileID);
                    // This file is specified using a UUID
                    // There should be an output connection to write it to
                    SocketChannel out = getOutputConnection(fileID);

                    if (out != null) {
                        try {
                            while (buffer.hasRemaining()) {
                                out.write(buffer);
                            }
                        } catch (IOException e) {
                            logger.warn("AMUSE output connection lost", e);
                            removeOutputConnection(fileID);
                        }
                    }
                }
            } catch (IOException e) {
                logger.error("Error while handling output", e);
            }

        }
    }

}
