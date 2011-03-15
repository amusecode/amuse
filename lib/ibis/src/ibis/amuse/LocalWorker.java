package ibis.amuse;

import ibis.deploy.Job;
import ibis.ipl.Ibis;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;
import ibis.util.ThreadPool;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;
import java.util.UUID;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Representation of worker at local machine running amuse script. Starts and connects to a "RemoteWorker" to perform the actual work.
 * 
 * @author Niels Drost
 * 
 */
public class LocalWorker implements Runnable {

    public static final String MAGIC_STRING = "magic_string";

    private static final Logger logger = LoggerFactory
            .getLogger(LocalWorker.class);

    private final SocketChannel channel;

    private final String codeName;

    private final String hostname;

    private final UUID id;

    private final Job job;

    private final ReceivePort receivePort;

    private final SendPort sendPort;

    /*
     * Initializes worker by reading settings from amuse, deploying the worker
     * process on a (possibly remote) machine, and waiting for a connection from
     * the worker
     */

    LocalWorker(SocketChannel socket, Ibis ibis, Deployment deployment)
            throws IOException {
        this.channel = socket;

        try {
            id = UUID.randomUUID();

            logger.info("New connection from "
                    + socket.socket().getRemoteSocketAddress());

            // read magic string, to make sure we are talking to amuse

            ByteBuffer magic = ByteBuffer.allocate(MAGIC_STRING.length());

            while (magic.hasRemaining()) {
                long read = channel.read(magic);

                if (read == -1) {
                    throw new IOException(
                            "Connection closed on reading magic string");
                }
            }

            String receivedString = new String(magic.array(), "UTF-8");
            if (!receivedString.equals(MAGIC_STRING)) {
                throw new IOException("magic string " + MAGIC_STRING
                        + " not received. Instead got: " + receivedString);
            }

            // read initialization call

            AmuseMessage initRequest = new AmuseMessage();
            initRequest.readFrom(channel);

            if (initRequest.getFunctionID() != AmuseMessage.FUNCTION_ID_INIT) {
                throw new IOException(
                        "first call to worker must be init function");
            }

            codeName = initRequest.getString(0);
            hostname = initRequest.getString(1);

            // initialize ibis ports

            receivePort = ibis
                    .createReceivePort(Daemon.portType, id.toString());
            receivePort.enableConnections();

            sendPort = ibis.createSendPort(Daemon.portType);

            // start deployment of worker (possibly on remote machine)

            job = deployment.deploy(codeName, hostname, id);

            // we expect a "hello" message from the worker.
            ReadMessage readMessage = receivePort.receive();
            ReceivePortIdentifier remotePort = (ReceivePortIdentifier) readMessage
                    .readObject();
            readMessage.finish();

            // connect to the worker
            sendPort.connect(remotePort);

            // send OK reply to amuse

            AmuseMessage initReply = new AmuseMessage(initRequest.getCallID(),
                    initRequest.getFunctionID(), initRequest.getCount());

            initReply.writeTo(channel);

            ThreadPool.createNew(this, "Amuse worker");

        } catch (Exception e) {
            // report error to amuse
            AmuseMessage errormessage = new AmuseMessage(0,
                    AmuseMessage.FUNCTION_ID_INIT, 1, e.getMessage());
            errormessage.writeTo(channel);
            throw new IOException("error initializing code", e);
        }
    }

    public void run() {
        boolean running = true;

        while (running) {
            AmuseMessage request = new AmuseMessage();
            AmuseMessage result = new AmuseMessage();

            try {
                request.readFrom(channel);

                if (request.getFunctionID() == AmuseMessage.FUNCTION_ID_STOP) {
                    // this will be the last call we perform
                    running = false;
                }

                if (request.getFunctionID() == AmuseMessage.FUNCTION_ID_REDIRECT_OUTPUT) {
                    logger.warn("Redirect output function not supported by IbisChannel");
                }

                WriteMessage writeMessage = sendPort.newMessage();
                request.writeTo(writeMessage);
                writeMessage.flush();

                ReadMessage readMessage = receivePort.receive();
                result.readFrom(readMessage);
                readMessage.finish();

                if (result.getError() != null) {
                    logger.warn("Error while doing call at worker",
                            result.getError());
                }

                // forward result to the channel
                result.writeTo(channel);

            } catch (IOException e) {
                logger.error("Error on handling call", e);
                // report error to amuse
                AmuseMessage errormessage = new AmuseMessage(request.getCallID(),
                        request.getFunctionID(), request.getCount(), e.getMessage());
                try {
                    errormessage.writeTo(channel);
                } catch (IOException e1) {
                    logger.error(
                            "Error while returning error message to amuse", e1);
                }
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e1) {
                    //IGNORE
                }
            }
        }
        end();
        logger.info(this + " done!");
    }

    private void end() {
        job.kill();
        
        try {
            sendPort.close();
        } catch (IOException e) {
            logger.error("Error closing sendport", e);
        }
        
        try {
            receivePort.close(1000);
        } catch (IOException e) {
            logger.error("Error closing receiveport", e);
        }
    }

    public String toString() {
        return "Worker \"" + id + "\" running \"" + codeName + "@" + hostname + "\"";
    }
}
