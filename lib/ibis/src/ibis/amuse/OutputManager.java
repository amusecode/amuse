package ibis.amuse;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.ipl.Ibis;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;

/**
 * Takes care of getting all output to its final destination. Either a local
 * file, or the AMUSE process.
 * 
 * @author Niels Drost
 * 
 */
public class OutputManager extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(OutputManager.class);

    public static final int MAX_MESSAGE_SIZE = 10240;

    private final Ibis ibis;
    private final ReceivePort receivePort;

    private final Map<UUID, SocketChannel> outputConnections;

    OutputManager(Ibis ibis) throws IOException {
        this.ibis = ibis;

        outputConnections = new HashMap<UUID, SocketChannel>();

        this.receivePort = ibis.createReceivePort(Daemon.portType, "output");
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
                    logger.debug("Got " + count + " bytes for file " + file);
                    // not a UUID, just write it to a file
                    FileOutputStream out = new FileOutputStream(file, true);
                    out.getChannel().write(buffer);
                    out.flush();
                    out.close();
                } else {
                    logger.debug("Got " + count + " bytes for stream " + fileID);
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
