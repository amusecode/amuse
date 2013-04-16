package ibis.amuse;

/**
 * Class representing a code that is used via JNI.
 * 
 * @author Niels Drost
 * 
 */
public class JNICodeInterface implements CodeInterface {

    private final WorkerInfo info;

    JNICodeInterface(WorkerInfo info) throws Exception {
        this.info = info;
    }

    public void init() throws CodeException {
        try {
            String library = info.getCodeName();

            System.loadLibrary(library);

            nativeInit(info.getCodeName());
        } catch (Throwable t) {
            throw new CodeException("Could not load worker library", t);
        }
    }

    /**
     * Sets the message used as input in a call. Also used when the buffers
     * inside the message change
     * 
     * @param message
     *            the message
     */
    public native void setRequestMessage(AmuseMessage message);

    /**
     * Sets the message used as output in a call
     * 
     * @param message
     *            the message
     */
    public native void setResultMessage(AmuseMessage message);

    /**
     * Performs a call. Uses data from the request and result message
     * 
     * @throws CodeException
     *             if the function does not exist, or any other error occurs.
     */
    public native void call() throws CodeException;

    /**
     * Initialize native code, if needed
     * 
     * @throws Exception
     *             if initialization fails
     */
    private native void nativeInit(String codeName) throws Exception;

    @Override
    public void end() {
    }

    @Override
    public int getResult() {
        return 0;
    }

}
