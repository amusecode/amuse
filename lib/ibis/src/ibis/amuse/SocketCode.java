package ibis.amuse;

import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;

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

    private final File executable;

    // messages shared with native code to perform actual calls

    private final AmuseMessage requestMessage;
    private final AmuseMessage resultMessage;

    // Ibis communication stuff

    private final ReceivePort receivePort;
    private final SendPort sendPort;

    SocketCode(String codeName, String codeDir, ReceivePort receivePort, SendPort sendPort)
            throws IOException {
        this.receivePort = receivePort;
        this.sendPort = sendPort;

        requestMessage = new AmuseMessage();
        resultMessage = new AmuseMessage();
        
        executable = new File(codeDir + File.separator + codeName);

        if (!executable.canExecute()) {
            throw new IOException("Cannot find executable for code " + codeName + ": " + executable);
        }
    }

    @Override
    /**
     * Continuously receives a message, performs a call, sends a reply.
     */
    public void run() {
        boolean running = true;

        while (running) {
            try {
                logger.debug("Receiving call message");
                ReadMessage readMessage = receivePort.receive();
                
                logger.debug("Reading call request from IPL message");

                boolean changed = requestMessage.readFrom(readMessage);
                
                readMessage.finish();
                
                               int functionID = requestMessage.getFunctionID();
                
                if (functionID == AmuseMessage.FUNCTION_ID_STOP) {
                    //final request handled
                    running = false;
                }

               
                
                logger.debug("Performing call for function " + functionID);
                try {
                    // perform call. Will put result in result message
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
            } catch (Throwable e) {
                logger.error("Error while handling request", e);
                try {
                Thread.sleep(1000);
                } catch (Exception e2) {
                    //IGNORE
                }
            }
        }

    }
    
//    
//    public static void main(String[] arguments) throws Exception {
//        SocketCode code = new SocketCode(arguments[0], null, null);
//        
//        code.requestMessage.clear();
//        code.requestMessage.setFunctionID(0);
//        code.requestMessage.setCallCount(1);
//
//        
//    }
}
