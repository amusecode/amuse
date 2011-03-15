package ibis.amuse;

import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class representing a code. Also provides the (native) interface to it.
 * 
 * @author Niels Drost
 * 
 */
public class CommunityCode implements Runnable {

    private static final Logger logger = LoggerFactory
            .getLogger(CommunityCode.class);

    private final String codeName;

    // messages shared with native code to perform actual calls

    private final AmuseMessage requestMessage;
    private final AmuseMessage resultMessage;

    // Ibis communication stuff

    private final ReceivePort receivePort;
    private final SendPort sendPort;

    // native functions

    /**
     * Returns the maximum result size of a single call per type
     */
    private native int[] getMaxResultSize();

    /**
     * Performs a call. Uses data from the request and result message
     * 
     * @throws Exception
     *             if the function does not exist, or any other error occurs.
     */
    private native void call() throws Exception;

    /**
     * Sets the message used as input in a call. Also used when the buffers
     * inside the message change
     * 
     * @param message
     *            the message
     */
    private native void setRequestMessage(AmuseMessage message);

    /**
     * Sets the message used as output in a call
     * 
     * @param message
     *            the message
     */
    private native void setResultMessage(AmuseMessage message);

    CommunityCode(String codeName, ReceivePort receivePort, SendPort sendPort) {
        this.codeName = codeName;
        this.receivePort = receivePort;
        this.sendPort = sendPort;

        System.loadLibrary("ibis-amuse-" + codeName);

        requestMessage = new AmuseMessage();
        resultMessage = new AmuseMessage();

        setRequestMessage(requestMessage);
        setResultMessage(resultMessage);
    }

    @Override
    /**
     * Continuously receives a message, performs a call, sends a reply.
     */
    public void run() {
        boolean running = true;
        
        while (running) {
            try {
                ReadMessage readMessage = receivePort.receive();

                if (requestMessage.readFrom(readMessage)) {
                    // read method indicates one or more buffers has changed in
                    // the message
                    setRequestMessage(requestMessage);
                }
                
                readMessage.finish();
                
                if (requestMessage.getFunctionID() == AmuseMessage.FUNCTION_ID_STOP) {
                    running = false;
                }
                
                try {
                    //perform call. Will put result in result message
                    call();
                } catch (Exception e) {
                    //put an exception in the result message
                    resultMessage.clear();
                    resultMessage.setCallID(requestMessage.getCallID());
                    resultMessage.setFunctionID(requestMessage.getFunctionID());
                    resultMessage.setCount(requestMessage.getCount());
                    resultMessage.setError(e.getMessage());
                    
                }
                
                WriteMessage writeMessage = sendPort.newMessage();
                
                resultMessage.writeTo(writeMessage);
                
                writeMessage.finish();
            } catch (IOException e) {
                logger.error("Error while handling request", e);
            }
        }

    }
}
