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
import ibis.ipl.IbisIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Forwards output of a worker to the "master" Ibis process.
 * 
 * @author Niels Drost
 *
 */
public class OutputForwarder extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(OutputForwarder.class);

    private final InputStream input;

    private final String file;

    private final Ibis ibis;

    public OutputForwarder(InputStream input, String file, Ibis ibis) throws IOException {
        this.input = new BufferedInputStream(input);
        this.file = file;
        this.ibis = ibis;

        setDaemon(false);
        setName("output forwarder");
        start();
    }

    /**
     * Forwards input stream to given output stream.
     */
    public void run() {
        byte[] buffer = new byte[WorkerOutputManager.MAX_MESSAGE_SIZE];
        SendPort sendPort = null;

        logger.debug("Forwarding output to " + file);

        try {

            if (file.equalsIgnoreCase("/dev/null")) {
                sendPort = null;
                logger.info("Discarding code output");
            } else {
                logger.debug("Creating sendport");
                sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);
                logger.debug("getting daemon IbisIdentifier");
                IbisIdentifier daemon = ibis.registry().getElectionResult("amuse");
                logger.debug("Connecting to output manager port");
                sendPort.connect(daemon, "output");
            }

            logger.debug("Starting with forwarding output");

            while (true) {
                int count = input.read(buffer);

                if (count == -1) {
                    // we're done (final block will close sendport)
                    return;
                }

                if (sendPort != null) {
                    logger.trace("Sending message with {} bytes", count);
                    WriteMessage message = sendPort.newMessage();
                    message.writeString(file);
                    message.writeInt(count);
                    message.writeArray(buffer, 0, count);
                    message.finish();
                }
            }
        } catch (IOException error) {
            if (!error.getMessage().equals("Stream closed")) {
                logger.error("Error on forwarding code output", error);
            }
        } finally {
            if (sendPort != null) {
                try {
                    sendPort.close();
                } catch (IOException e) {
                    // IGNORE
                }
            }
        }

    }

}
