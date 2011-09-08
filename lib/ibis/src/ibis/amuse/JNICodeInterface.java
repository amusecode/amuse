package ibis.amuse;

import java.io.IOException;

/**
 * Class representing a code that is used via JNI.
 * 
 * @author Niels Drost
 * 
 */
public class JNICodeInterface extends CodeInterface {

    JNICodeInterface(String id, PoolInfo poolInfo, String codeName)
            throws Exception {
        super(id, poolInfo);

        AmuseMessage initRequest = receiveInitRequest();

        try {
            String library = codeName;

            System.loadLibrary(library);

            init(codeName);
        } catch (Throwable t) {
            IOException error = new IOException(
                    "Could not load worker library", t);
            sendInitReply(initRequest.getCallID(), error);
            
            end();
                
            throw error;
        }
        sendInitReply(initRequest.getCallID());
    }

    /**
     * Sets the message used as input in a call. Also used when the buffers
     * inside the message change
     * 
     * @param message
     *            the message
     */
    native void setRequestMessage(AmuseMessage message);

    /**
     * Sets the message used as output in a call
     * 
     * @param message
     *            the message
     */
    native void setResultMessage(AmuseMessage message);

    /**
     * Performs a call. Uses data from the request and result message
     * 
     * @throws Exception
     *             if the function does not exist, or any other error occurs.
     */
    native void call() throws Exception;

    /**
     * Initialize native code, if needed
     * 
     * @throws Exception
     *             if initialization fails
     */
    private native void init(String codeName) throws Exception;

}
