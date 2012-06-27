package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class OutputForwarder extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(OutputForwarder.class);

    private final InputStream input;

    private final String file;

    private final Ibis ibis;

    /**
     * @param input
     *            Input stream
     * @param output
     *            Stream to forward output to
     * @param codeName
     *            Prefix to add to all lines of output
     * 
     * @throws IOException
     *             if the reading stream cannot be created.
     */
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
        byte[] buffer = new byte[OutputManager.MAX_MESSAGE_SIZE];
        SendPort sendPort = null;

        logger.debug("Forwarding output to " + file);

        try {

            if (file.equalsIgnoreCase("/dev/null")) {
                sendPort = null;
                logger.info("Discarding code output");
            } else {
                logger.debug("Creating sendport");
                sendPort = ibis.createSendPort(Daemon.portType);
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
                    logger.debug("Sending message with " + count + " bytes");
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
