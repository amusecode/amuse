package ibis.amuse;

interface CodeInterface {

    /**
     * Initialize.
     * 
     * @throws CodeException
     *             if initialization fails
     */
    void init() throws CodeException;

    /**
     * Sets the message used as input in a call. Also used when the buffers
     * inside the message change
     * 
     * @param message
     *            the message
     */
    void setRequestMessage(AmuseMessage message) throws CodeException;

    /**
     * Sets the message used as output in a call
     * 
     * @param message
     *            the message
     */
    void setResultMessage(AmuseMessage message) throws CodeException;

    /**
     * Performs a call. Uses data from the request and result message
     * 
     * @throws CodeException
     *             if the function does not exist, or any other error occurs.
     */
    void call() throws CodeException;

    void end();

    int getResult();

}
