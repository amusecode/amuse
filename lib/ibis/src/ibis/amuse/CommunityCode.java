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

    /**
     * Performs a call. Uses data from the request and result message
     * 
     * @throws Exception
     *             if the function does not exist, or any other error occurs.
     */
    private native void call() throws Exception;

    /**
     * Initialize native code, if needed
     * 
     * @throws Exception if initialization fails
     */
    private native void init(String codeName) throws Exception;
    
    CommunityCode(String codeName, ReceivePort receivePort, SendPort sendPort)
            throws IOException {
        this.codeName = codeName;
        this.receivePort = receivePort;
        this.sendPort = sendPort;

        requestMessage = new AmuseMessage();
        resultMessage = new AmuseMessage();

        String library = codeName;
        
        //this is a MPI worker, load mpi worker lib
        if (codeName.contains("/")) {
        	logger.info("MPI worker mode...");
        	library = "mpi_ibis_worker";
        	
//            System.loadLibrary("mpi");

        }

        try {
            System.loadLibrary(library);

            init(codeName);
            setRequestMessage(requestMessage);
            setResultMessage(resultMessage);
        } catch (Throwable t) {
            logger.error("Could not load worker library \"" + library + "\"", t);
            throw new IOException("Could not load worker library", t);
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

                if (requestMessage.readFrom(readMessage)) {
                    // read method indicates one or more buffers has changed in
                    // the message
                    logger.debug("(re)setting request message");
                    setRequestMessage(requestMessage);
                }

                readMessage.finish();

                int functionID = requestMessage.getFunctionID();
                
                if (functionID == AmuseMessage.FUNCTION_ID_STOP) {
                    //final request handled
                    running = false;
                }

                //clean message
                resultMessage.clear();
                
                logger.debug("Performing call for function " + functionID);
                try {
                    // perform call. Will put result in result message
                    call();
                } catch (Exception e) {
                    logger.error("exception on performing call", e);
                    // put an exception in the result message
                    resultMessage.clear();
                    resultMessage.setCallID(requestMessage.getCallID());
                    resultMessage.setFunctionID(requestMessage.getFunctionID());
                    resultMessage.setCount(requestMessage.getCount());
                    resultMessage.setError(e.getMessage());
                }
                
                logger.debug("result: " + resultMessage);

                WriteMessage writeMessage = sendPort.newMessage();

                resultMessage.writeTo(writeMessage);

                writeMessage.finish();
                
                logger.debug("Done performing call for function " + functionID);
            } catch (IOException e) {
                logger.error("Error while handling request", e);
            }
        }

    }
}
