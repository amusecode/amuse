package ibis.amuse;

import ibis.ipl.ReadMessage;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.io.Serializable;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AmuseMessage implements Serializable {

    public static final int HEADER_SIZE = 10; // integers

    // 4 byte flags field.
    public static final int HEADER_FLAGS = 0;

    // content of flags field (first 4 bytes of message header) currently:
    // - endianness
    // - if an exception has occurred
    public static final int ENDIANNES_FLAG = 0;
    public static final int ERROR_FLAG = 1;

    public static final int HEADER_CALL_ID_INDEX = 1;
    public static final int HEADER_FUNCTION_ID_INDEX = 2;
    public static final int HEADER_CALL_COUNT_INDEX = 3;
    public static final int HEADER_INT_COUNT_INDEX = 4;
    public static final int HEADER_LONG_COUNT_INDEX = 5;
    public static final int HEADER_FLOAT_COUNT_INDEX = 6;
    public static final int HEADER_DOUBLE_COUNT_INDEX = 7;
    public static final int HEADER_BOOLEAN_COUNT_INDEX = 8;
    public static final int HEADER_STRING_COUNT_INDEX = 9;

    public static final int SIZEOF_INT = 4;
    public static final int SIZEOF_LONG = 8;
    public static final int SIZEOF_FLOAT = 8;
    public static final int SIZEOF_DOUBLE = 4;
    public static final int SIZEOF_BOOLEAN = 1;

    public static final byte TRUE_BYTE = (1 & 0xFF);
    public static final byte FALSE_BYTE = (0 & 0xFF);

    // "random" numbers denoting big/little endian
    public static final byte BIG_ENDIAN = 67;
    public static final byte LITTLE_ENDIAN = 101;

    public static final int FUNCTION_ID_INIT = 10101010;
    public static final int FUNCTION_ID_STOP = 0;
    public static final int FUNCTION_ID_REDIRECT_OUTPUT = 1141573512;

    private static final long serialVersionUID = 1L;

    private static final Logger logger = LoggerFactory
            .getLogger(AmuseMessage.class);

    private static boolean hasRemaining(ByteBuffer... buffers) {
        for (ByteBuffer buffer : buffers) {
            if (buffer.hasRemaining()) {
                return true;
            }
        }
        return false;
    }

    static void readAll(SocketChannel channel, ByteBuffer... bytes)
            throws IOException {

        while (hasRemaining(bytes)) {
            long read = channel.read(bytes);

            if (read == -1) {
                throw new IOException("Connection closed on reading data");
            }
        }
    }

    private final ByteBuffer headerBytes;

    private ByteBuffer intBytes;

    private ByteBuffer longBytes;

    private ByteBuffer floatBytes;

    private ByteBuffer doubleBytes;

    private ByteBuffer booleanBytes;

    private ByteBuffer stringHeaderBytes;

    private ByteBuffer[] byteBuffers;
    private ByteBuffer[] allButStringByteBuffers;

    // UTF-8 encoded strings
    private ByteBuffer[] stringBytes;

    // view of buffers (for easy access)

    private final IntBuffer header;

    private IntBuffer stringHeader;

    /**
     * Empty message.
     */
    public AmuseMessage() {
        headerBytes = ByteBuffer.allocateDirect(HEADER_SIZE * SIZEOF_INT);
        intBytes = ByteBuffer.allocateDirect(0);
        longBytes = ByteBuffer.allocateDirect(0);
        floatBytes = ByteBuffer.allocateDirect(0);
        doubleBytes = ByteBuffer.allocateDirect(0);
        booleanBytes = ByteBuffer.allocateDirect(0);
        stringHeaderBytes = ByteBuffer.allocateDirect(0);
        stringBytes = new ByteBuffer[0];

        allButStringByteBuffers = new ByteBuffer[] { headerBytes, intBytes,
                longBytes, floatBytes, doubleBytes, booleanBytes,
                stringHeaderBytes };

        // no string buffers yet
        byteBuffers = allButStringByteBuffers;

        header = headerBytes.asIntBuffer();
        stringHeader = stringHeaderBytes.asIntBuffer();
    }

    AmuseMessage(int callID, int functionID, int count) {
        this();

        setCallID(callID);
        setFunctionID(functionID);
        setCount(count);
    }

    /**
     * Massages with an exception
     * 
     * @param callID
     *            id of the call that generated the exception
     * @param functionID
     *            function id tried
     * @param error
     *            a description of the error that occurred
     */
    AmuseMessage(int callID, int functionID, int count, String error) {
        this();

        setCallID(callID);
        setFunctionID(functionID);
        setCount(count);
        setError(error);
    }

    public void clear() {
        headerBytes.clear();
        // stuff full of zeros
        headerBytes.put(new byte[headerBytes.capacity()]);
        setByteOrder(ByteOrder.nativeOrder());
    }

    /**
     * Change the byte order of this message.
     * 
     * @param order
     *            The new byte-order
     */
    private void setByteOrder(ByteOrder order) {
        if (order == ByteOrder.BIG_ENDIAN) {
            headerBytes.put(BIG_ENDIAN);
        } else {
            headerBytes.put(LITTLE_ENDIAN);
        }

        for (ByteBuffer buffer : getByteBuffers(false)) {
            buffer.order(order);
        }
    }

    /**
     * Change the byte order of this message. Also swaps the content of all the
     * buffers, if requested.
     * 
     * @param order
     *            The new byte-order
     * @param swapContent
     *            if True, all data contained in this message is byte-order
     *            swapped.
     * @throws IOException
     *             if the byte order cannot be determined
     */
    void setByteOrder(ByteOrder order, boolean swapContent) throws IOException {
        ByteOrder oldOrder = getByteOrder();

        if (order == oldOrder) {
            // done! :-)
            return;
        }

        throw new IOException("byte swapping not implemented yet!");
    }

    private ByteOrder getByteOrder() {
        if (headerBytes.get(0) == BIG_ENDIAN) {
            return ByteOrder.BIG_ENDIAN;
        } else if (headerBytes.get(0) == LITTLE_ENDIAN) {
            return ByteOrder.LITTLE_ENDIAN;
        } else {
            throw new RuntimeException("endiannes not specified in header");
        }
    }

    public void setCount(int count) {
        header.put(HEADER_CALL_COUNT_INDEX, count);
    }

    public void setFunctionID(int functionID) {
        header.put(HEADER_FUNCTION_ID_INDEX, functionID);
    }

    public void setCallID(int callID) {
        header.put(HEADER_CALL_ID_INDEX, callID);
    }

    public void setError(String error) {
        // clear data from message
        header.put(HEADER_INT_COUNT_INDEX, 0);
        header.put(HEADER_LONG_COUNT_INDEX, 0);
        header.put(HEADER_FLOAT_COUNT_INDEX, 0);
        header.put(HEADER_DOUBLE_COUNT_INDEX, 0);
        header.put(HEADER_BOOLEAN_COUNT_INDEX, 0);
        header.put(HEADER_STRING_COUNT_INDEX, 1);

        // set error state
        headerBytes.put(ERROR_FLAG, TRUE_BYTE);

        ensurePrimitiveCapacity();

        try {
            // set first string to exception message
            byte[] bytes;
            
            bytes = error.getBytes("UTF-8");

            stringHeader.put(0, bytes.length);

            ensureStringsCapacity();

            stringBytes[0].clear();
            stringBytes[0].put(bytes);
        } catch (UnsupportedEncodingException e) {
            logger.error("could not set error", e);
            stringHeader.put(0, 0);
        }
    }

    public boolean isErrorState() {
        return headerBytes.get(ERROR_FLAG) == TRUE_BYTE;
    }

    public String getError() throws IOException {
        if (!isErrorState()) {
            return null;
        }
        return getString(0);
    }

    public int getCallID() {
        return header.get(HEADER_CALL_ID_INDEX);
    }

    public int getFunctionID() {
        return header.get(HEADER_FUNCTION_ID_INDEX);
    }

    public int getCount() {
        return header.get(HEADER_CALL_COUNT_INDEX);
    }

    public int getIntCount() {
        return header.get(HEADER_STRING_COUNT_INDEX);
    }

    public int getLongCount() {
        return header.get(HEADER_LONG_COUNT_INDEX);
    }

    public int getFloatCount() {
        return header.get(HEADER_FLOAT_COUNT_INDEX);
    }

    public int getDoubleCount() {
        return header.get(HEADER_DOUBLE_COUNT_INDEX);
    }

    public int getBooleanCount() {
        return header.get(HEADER_BOOLEAN_COUNT_INDEX);
    }

    public int getStringCount() {
        return header.get(HEADER_STRING_COUNT_INDEX);
    }

    public String getString(int index) throws IOException {
        if (getStringCount() <= index) {
            throw new IOException("cannot get string at index " + index
                    + " in call" + this);
        }

        if (stringBytes.length <= index) {
            throw new IOException("cannot get string at index " + index
                    + " in call" + this + " header does not match content!");

        }

        int utf8length = stringHeader.get(index);

        if (stringBytes[index].hasArray()) {
            return new String(stringBytes[index].array(), 0, utf8length,
                    "UTF-8");
        }
        byte[] bytes = new byte[utf8length];
        stringBytes[index].position(0);
        stringBytes[index].limit(utf8length);
        stringBytes[index].get(bytes);

        return new String(stringBytes[index].array(), 0, utf8length, "UTF-8");
    }

    /**
     * Get all buffers, possibly including the buffers containing the strings.
     * 
     * @return all buffers.
     * 
     * @param includeStringBuffers
     *            if true, the buffers for holding the values of strings will be
     *            included.
     */
    public ByteBuffer[] getByteBuffers(boolean includeStringBuffers) {
        if (includeStringBuffers) {
            return byteBuffers;
        } else {
            return allButStringByteBuffers;
        }
    }

    private void setPrimitiveLimitsFromHeader() throws IOException {
        intBytes.clear().limit(getIntCount() * SIZEOF_INT);
        longBytes.clear().limit(getLongCount() * SIZEOF_LONG);
        floatBytes.clear().limit(getFloatCount() * SIZEOF_FLOAT);
        doubleBytes.clear().limit(getDoubleCount() * SIZEOF_DOUBLE);
        booleanBytes.clear().limit(getBooleanCount() * SIZEOF_BOOLEAN);
        stringHeaderBytes.clear().limit(getStringCount() * SIZEOF_INT);
    }

    private void setStringLimitsFromHeader() throws IOException {
        if (getStringCount() > stringBytes.length) {
            throw new IOException(
                    "Amuse message in inconsistent state, strign count greater than number of string buffers");
        }

        for (int i = 0; i < getStringCount(); i++) {
            int utf8Length = stringHeader.get(i);

            stringBytes[i].clear().limit(utf8Length);
        }

        // set the limit of the rest of the string bytes to 0
        for (int i = getStringCount(); i < stringBytes.length; i++) {
            stringBytes[i].limit(0);
        }
    }

    void writeTo(SocketChannel channel) throws IOException {
        headerBytes.clear();
        setPrimitiveLimitsFromHeader();
        setStringLimitsFromHeader();

        // write to channel
        channel.write(byteBuffers);
    }

    void writeTo(WriteMessage writeMessage) throws IOException {
        headerBytes.clear();
        setPrimitiveLimitsFromHeader();
        setStringLimitsFromHeader();

        for (ByteBuffer buffer : byteBuffers) {
            writeMessage.writeByteBuffer(buffer);
        }
    }

    // make sure there is enough space for each primitive buffer
    // (including the string header)
    private boolean ensurePrimitiveCapacity() {
        boolean buffersUpdated = false;

        if (getIntCount() < (intBytes.capacity() * SIZEOF_INT)) {
            intBytes = ByteBuffer.allocateDirect(getIntCount() * SIZEOF_INT);
            intBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getLongCount() < (longBytes.capacity() * SIZEOF_LONG)) {
            longBytes = ByteBuffer.allocateDirect(getLongCount() * SIZEOF_LONG);
            longBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getFloatCount() < (floatBytes.capacity() * SIZEOF_FLOAT)) {
            floatBytes = ByteBuffer.allocateDirect(getFloatCount()
                    * SIZEOF_FLOAT);
            floatBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getDoubleCount() < (doubleBytes.capacity() * SIZEOF_DOUBLE)) {
            doubleBytes = ByteBuffer.allocateDirect(getDoubleCount()
                    * SIZEOF_DOUBLE);
            doubleBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getBooleanCount() < (booleanBytes.capacity() * SIZEOF_BOOLEAN)) {
            booleanBytes = ByteBuffer.allocateDirect(getBooleanCount()
                    * SIZEOF_BOOLEAN);
            booleanBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getStringCount() < (stringHeaderBytes.capacity() * SIZEOF_INT)) {
            stringHeaderBytes = ByteBuffer.allocateDirect(getStringCount()
                    * SIZEOF_INT);
            stringHeaderBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (buffersUpdated) {
            allButStringByteBuffers = new ByteBuffer[] { headerBytes, intBytes,
                    longBytes, floatBytes, doubleBytes, booleanBytes,
                    stringHeaderBytes };

            if (logger.isDebugEnabled()) {
                logger.debug("Updated buffers to " + this);
            }
        }

        return buffersUpdated;
    }

    private boolean ensureStringsCapacity() {
        // checking if the string header is big enough is checked above, so we
        // only check if all strings listed in the header
        boolean buffersUpdated = false;

        if (stringBytes.length < getStringCount()) {
            ByteBuffer[] oldStringBytes = stringBytes;
            stringBytes = new ByteBuffer[getStringCount()];
            for (int i = 0; i < oldStringBytes.length; i++) {
                stringBytes[i] = oldStringBytes[i];
            }
            buffersUpdated = true;
        }

        for (int i = 0; i < getStringCount(); i++) {
            if (stringBytes[i] == null
                    || stringBytes[i].capacity() < stringHeader.get(i)) {
                stringBytes[i] = ByteBuffer.allocateDirect(stringHeader.get(i));
                buffersUpdated = true;
            }
        }
        return buffersUpdated;
    }

    boolean readFrom(SocketChannel channel) throws IOException {
        boolean updatedBuffers = false;

        logger.debug("receiving header from channel");

        headerBytes.clear();

        readAll(channel, headerBytes);

        // set buffers to byte order specified in buffer
        setByteOrder(getByteOrder());

        if (ensurePrimitiveCapacity()) {
            updatedBuffers = true;
        }

        // then, set limits for primitive buffers, and receive those

        setPrimitiveLimitsFromHeader();

        logger.debug("receiving primitives from channel");

        // we also request to read the header, but its position is already
        // equal to its limit, so no bytes are read into it.
        readAll(channel, allButStringByteBuffers);

        // make sure there is enough space for the strings
        if (ensureStringsCapacity()) {
            updatedBuffers = true;
        }

        // set the limits
        setStringLimitsFromHeader();

        logger.debug("receiving strings from channel");

        // and receive!
        readAll(channel, stringBytes);

        if (logger.isDebugEnabled()) {
            logger.debug("done receiving message from channel: " + this);
        }

        return updatedBuffers;
    }

    public boolean readFrom(ReadMessage readMessage) throws IOException {
        boolean updatedBuffers = false;

        logger.debug("reading header from message");

        headerBytes.clear();

        readMessage.readByteBuffer(headerBytes);

        // set buffers to byte order specified in buffer
        setByteOrder(getByteOrder());

        if (ensurePrimitiveCapacity()) {
            updatedBuffers = true;
        }

        // then, set limits for primitive buffers, and receive those

        setPrimitiveLimitsFromHeader();

        logger.debug("reading primitives from message");

        readMessage.readByteBuffer(intBytes);
        readMessage.readByteBuffer(longBytes);
        readMessage.readByteBuffer(floatBytes);
        readMessage.readByteBuffer(doubleBytes);
        readMessage.readByteBuffer(booleanBytes);
        readMessage.readByteBuffer(stringHeaderBytes);

        // make sure there is enough space for the strings
        if (ensureStringsCapacity()) {
            updatedBuffers = true;
        }

        // set the limits
        setStringLimitsFromHeader();

        logger.debug("reading strings from message");

        // and receive!
        for (ByteBuffer buffer : stringBytes) {
            readMessage.readByteBuffer(buffer);
        }

        if (logger.isDebugEnabled()) {
            logger.debug("done receiving message from channel: " + this);
        }

        return updatedBuffers;
    }

    public String toString() {
        return "Call: function ID:" + getFunctionID() + " count:" + getCount()
                + " ints:" + getIntCount() + " longs: " + getLongCount()
                + " floats:" + getFloatCount() + " doubles:" + getDoubleCount()
                + " booleans:" + getBooleanCount() + " strings:"
                + getStringCount();
    }

}
