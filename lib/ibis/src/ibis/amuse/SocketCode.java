package ibis.amuse;

import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;

import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class representing a code that is used via a loopback socket interface.
 * 
 * @author Niels Drost
 * 
 */
public class SocketCode implements Runnable {

    private static final Logger logger = LoggerFactory
            .getLogger(SocketCode.class);

    private static final int ACCEPT_TRIES = 5;
    private static final int ACCEPT_TIMEOUT = 500; //ms

    private final File executable;

    // messages used for receiving/sending requests

    private final AmuseMessage requestMessage;
    private final AmuseMessage resultMessage;

    // Ibis communication stuff

    private final ReceivePort receivePort;
    private final SendPort sendPort;

    // local socket communication stuff

    private final ServerSocketChannel serverSocket;

    private final SocketChannel socket;

    private final Process process;

    SocketCode(String codeName, String codeDir, String amuseHome, ReceivePort receivePort,
            SendPort sendPort) throws IOException {
        this.receivePort = receivePort;
        this.sendPort = sendPort;

        requestMessage = new AmuseMessage();
        resultMessage = new AmuseMessage();

        executable = new File(amuseHome + File.separator + codeDir + File.separator + codeName);

        if (!executable.canExecute()) {
            throw new IOException("Cannot find executable for code " + codeName
                    + ": " + executable);
        }

        serverSocket = ServerSocketChannel.open();
        serverSocket.socket().bind(null);

        ProcessBuilder builder = new ProcessBuilder();
        builder.command(executable.toString(),
                Integer.toString(serverSocket.socket().getLocalPort()));

        process = builder.start();
        
        new OutputPrefixForwarder(process.getInputStream(), System.out, "stdout of " + codeName + ": ");
        new OutputPrefixForwarder(process.getErrorStream(), System.err, "stderr of " + codeName + ": ");

        logger.info("process started");

        socket = acceptConnection(serverSocket);
        
        logger.info("connection established");

    }
    
    private static SocketChannel acceptConnection(ServerSocketChannel serverSocket) throws IOException {
        serverSocket.configureBlocking(false);
        for(int i = 0; i < ACCEPT_TRIES;i++) {
            SocketChannel result = serverSocket.accept();
            
            if (result != null) {
                return result;
            }
            try {
                Thread.sleep(ACCEPT_TIMEOUT);
            } catch (InterruptedException e) {
                //IGNORE
            }
        }
        throw new IOException("worker not started, socket connection failed to initialize");
    }

    @Override
    /**
     * Continuously receives a message, performs a call, sends a reply.
     */
    public void run() {
        boolean running = true;
        long start, finish;

        while (running) {
            
            try {
                start = System.currentTimeMillis();
                logger.debug("Receiving call message");
                ReadMessage readMessage = receivePort.receive();

                logger.debug("Reading call request from IPL message");

                requestMessage.readFrom(readMessage);

                readMessage.finish();

                int functionID = requestMessage.getFunctionID();

                if (functionID == AmuseMessage.FUNCTION_ID_STOP) {
                    // final request handled
                    running = false;
                }

                logger.debug("Performing call for function " + functionID);
                try {
                    requestMessage.writeTo(socket);
                    resultMessage.readFrom(socket);
                } catch (Exception e) {
                    logger.error("exception on performing call", e);
                    // put an exception in the result message
                    resultMessage.clear();
                    resultMessage.setCallID(requestMessage.getCallID());
                    resultMessage.setFunctionID(requestMessage.getFunctionID());
                    resultMessage.setCallCount(requestMessage.getCount());
                    resultMessage.setError(e.getMessage());
                }

                logger.debug("result: " + resultMessage);

                WriteMessage writeMessage = sendPort.newMessage();

                resultMessage.writeTo(writeMessage);

                writeMessage.finish();

                logger.debug("Done performing call for function " + functionID);
                finish = System.currentTimeMillis();
                
                if (logger.isInfoEnabled()) {
                    logger.info("Call took " + (finish - start) + " ms");
                }
            } catch (Throwable e) {
                logger.error("Error while handling request", e);
                try {
                    Thread.sleep(1000);
                } catch (Exception e2) {
                    // IGNORE
                }
            }
        }
        process.destroy();
    }

    public static void main(String[] arguments) throws Exception {
        SocketCode code = new SocketCode(arguments[0], arguments[1], arguments[2], null, null);

        code.requestMessage.clear();
        code.requestMessage.setFunctionID(1644113439);
        code.requestMessage.setCallCount(1);

        code.requestMessage.writeTo(code.socket);
        code.resultMessage.readFrom(code.socket);

    }
}
