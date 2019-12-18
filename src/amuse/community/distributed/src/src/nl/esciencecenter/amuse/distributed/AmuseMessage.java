/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nl.esciencecenter.amuse.distributed;

import ibis.ipl.ReadMessage;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class AmuseMessage {

    private static final Logger logger = LoggerFactory.getLogger(AmuseMessage.class);

    public static final int HEADER_SIZE = 11; // integers

    // 4 byte flags field.
    public static final int HEADER_FLAGS = 0;

    // content of flags field (first 4 bytes of message header) currently:
    // - endianness
    // - if an exception has occurred
    public static final int HEADER_BIG_ENDIAN_FLAG = 0;
    public static final int HEADER_ERROR_FLAG = 1;

    public static final int HEADER_CALL_ID_INDEX = 1;
    public static final int HEADER_FUNCTION_ID_INDEX = 2;
    public static final int HEADER_CALL_COUNT_INDEX = 3;
    public static final int HEADER_INT_COUNT_INDEX = 4;
    public static final int HEADER_LONG_COUNT_INDEX = 5;
    public static final int HEADER_FLOAT_COUNT_INDEX = 6;
    public static final int HEADER_DOUBLE_COUNT_INDEX = 7;
    public static final int HEADER_BOOLEAN_COUNT_INDEX = 8;
    public static final int HEADER_STRING_COUNT_INDEX = 9;
    public static final int HEADER_UNITS_COUNT_INDEX = 10;

    public static final int SIZEOF_INT = 4;
    public static final int SIZEOF_LONG = 8;
    public static final int SIZEOF_FLOAT = 4;
    public static final int SIZEOF_DOUBLE = 8;
    public static final int SIZEOF_BOOLEAN = 1;

    public static final byte TRUE_BYTE = (1 & 0xFF);
    public static final byte FALSE_BYTE = (0 & 0xFF);

    public static final int FUNCTION_ID_INIT = 10101010;
    public static final int FUNCTION_ID_STOP = 0;
    public static final int FUNCTION_ID_REDIRECT_OUTPUT = 1141573512;

    private static boolean hasRemaining(ByteBuffer... buffers) {
        for (ByteBuffer buffer : buffers) {
            if (buffer.hasRemaining()) {
                return true;
            }
        }
        return false;
    }

    static void readAll(SocketChannel channel, ByteBuffer... bytes) throws IOException {

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

    private IntBuffer header;

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

        allButStringByteBuffers = new ByteBuffer[] { headerBytes, intBytes, longBytes, floatBytes, doubleBytes,
                booleanBytes, stringHeaderBytes };

        // no string buffers yet
        byteBuffers = allButStringByteBuffers;

        ByteOrder nativeOrder = ByteOrder.nativeOrder();

        for (ByteBuffer buffer : byteBuffers) {
            buffer.order(nativeOrder);
        }

        header = headerBytes.asIntBuffer();
        stringHeader = stringHeaderBytes.asIntBuffer();
    }

    public AmuseMessage(int callID, int functionID, int count) {
        this();

        setCallID(callID);
        setFunctionID(functionID);
        setCallCount(count);
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
    public AmuseMessage(int callID, int functionID, int count, String error) {
        this();

        setCallID(callID);
        setFunctionID(functionID);
        setCallCount(count);
        setError(error);
    }

    /**
     * @return the size of the data currently in this message (in bytes).
     */
    public long getDataSize() {
        long result = 0;

        result += getIntCount() * SIZEOF_INT;
        result += getLongCount() * SIZEOF_LONG;
        result += getFloatCount() * SIZEOF_FLOAT;
        result += getDoubleCount() * SIZEOF_DOUBLE;
        result += getBooleanCount() * SIZEOF_BOOLEAN;

        for (int i = 0; i < getStringCount(); i++) {
            result += stringHeader.get(i)+1 ; // account for zero
        }

        return result;
    }

    public void clear() {
        headerBytes.clear();

        // stuff full of zeros

        byte[] zeros = new byte[headerBytes.capacity()];

        // remember byte order
        zeros[HEADER_BIG_ENDIAN_FLAG] = headerBytes.get(HEADER_BIG_ENDIAN_FLAG);

        headerBytes.put(zeros);
    }

    /**
     * Change the byte order of this message.
     * 
     * @param order
     *            The new byte-order
     */
    private void setByteOrder(ByteOrder order) {
        if (order == ByteOrder.BIG_ENDIAN) {
            headerBytes.put(HEADER_BIG_ENDIAN_FLAG, TRUE_BYTE);
        } else {
            headerBytes.put(HEADER_BIG_ENDIAN_FLAG, FALSE_BYTE);
        }

        for (ByteBuffer buffer : getByteBuffers(false)) {
            buffer.order(order);
        }

        // re-create views, as the order-change may not become visible
        // otherwise
        headerBytes.clear();
        header = headerBytes.asIntBuffer();
        stringHeaderBytes.clear();
        stringHeader = stringHeaderBytes.asIntBuffer();
    }

    /**
     * Change the byte order of this message. Also swaps the content of all the buffers, if requested.
     * 
     * @param order
     *            The new byte-order
     * @param swapContent
     *            if True, all data contained in this message is byte-order swapped.
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
        if (headerBytes.get(HEADER_BIG_ENDIAN_FLAG) == TRUE_BYTE) {
            return ByteOrder.BIG_ENDIAN;
        } else if (headerBytes.get(HEADER_BIG_ENDIAN_FLAG) == FALSE_BYTE) {
            return ByteOrder.LITTLE_ENDIAN;
        } else {
            throw new RuntimeException("endiannes not specified in header");
        }
    }

    public void setCallCount(int count) {
        header.put(HEADER_CALL_COUNT_INDEX, count);
    }

    public void setFunctionID(int functionID) {
        header.put(HEADER_FUNCTION_ID_INDEX, functionID);
    }

    public void setCallID(int callID) {
        header.put(HEADER_CALL_ID_INDEX, callID);
    }

    public void setError(String error) {
        if (error == null) {
            error = "<empty>";
        }
        // clear data from message
        header.put(HEADER_INT_COUNT_INDEX, 0);
        header.put(HEADER_LONG_COUNT_INDEX, 0);
        header.put(HEADER_FLOAT_COUNT_INDEX, 0);
        header.put(HEADER_DOUBLE_COUNT_INDEX, 0);
        header.put(HEADER_BOOLEAN_COUNT_INDEX, 0);
        header.put(HEADER_STRING_COUNT_INDEX, 1);
        header.put(HEADER_UNITS_COUNT_INDEX, 0);

        // set error state
        headerBytes.put(HEADER_ERROR_FLAG, TRUE_BYTE);

        ensurePrimitiveCapacity();

        try {
            // set first string to exception message
            byte[] bytes;

            bytes = error.getBytes("UTF-8");

            stringHeader.put(0, bytes.length);

            ensureStringsCapacity(0);

            stringBytes[0].clear();
            stringBytes[0].put(bytes);
            stringBytes[0].put( (byte) 0); // add extra zero
        } catch (UnsupportedEncodingException e) {
            System.err.println("could not set error: " + e);
            stringHeader.put(0, 0);
        }
    }

    public void setBoolean(int index, boolean value) {
        if (value) {
            booleanBytes.put(index, TRUE_BYTE);
        } else {
            booleanBytes.put(index, FALSE_BYTE);
        }

    }

    public void addString(String value) {
        int position = header.get(HEADER_STRING_COUNT_INDEX);

        // add an extra string
        header.put(HEADER_STRING_COUNT_INDEX, position + 1);

        // make sure there is space in the header for the length of the
        // string
        ensurePrimitiveCapacity();

        // encode string to UTF-8
        byte[] bytes;

        try {
            if (value == null) {
                //set null values to an empty string
                bytes = new String().getBytes("UTF-8");
            } else {
                bytes = value.getBytes("UTF-8");
            }

            // set length of string in header
            stringHeader.put(position, bytes.length);

            // make sure there is space for the string
            ensureStringsCapacity(position);

            stringBytes[position].clear();
            stringBytes[position].put(bytes);
            stringBytes[position].put( (byte) 0); // add extra zero
        } catch (UnsupportedEncodingException e) {
            System.err.println("ERROR! UTF-8 not supported by the JVM!");
        }
    }

    public void setString(int index, String value) {
        // encode string to UTF-8
        byte[] bytes;

        try {
            if (value == null) {
                //set null values to an empty string
                bytes = new String().getBytes("UTF-8");
            } else {
                bytes = value.getBytes("UTF-8");
            }
            
            // set length of string in header
            stringHeader.put(index, bytes.length);

            // make sure there is space for the string
            ensureStringsCapacity(index);

            stringBytes[index].clear();
            stringBytes[index].put(bytes);
            stringBytes[index].put( (byte) 0); // add extra zero
        } catch (UnsupportedEncodingException e) {
            System.err.println("ERROR! UTF-8 not supported by the JVM!");
        }
    }

    public boolean isErrorState() {
        return headerBytes.get(HEADER_ERROR_FLAG) == TRUE_BYTE;
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

    public int getCallCount() {
        return header.get(HEADER_CALL_COUNT_INDEX);
    }

    public int getIntCount() {
        return header.get(HEADER_INT_COUNT_INDEX);
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

    public String getString(int index, boolean convertEmptyToNull) throws IOException {
        String result = getString(index);

        if (convertEmptyToNull && result.isEmpty()) {
            result = null;
        }
        return result;
    }

    public String getString(int index) throws IOException {
        if (getStringCount() <= index) {
            throw new IOException("cannot get string at index " + index + " in call" + this);
        }

        if (stringBytes.length <= index) {
            throw new IOException("cannot get string at index " + index + " in call" + this 
                     + " header does not match content!");

        }

        int utf8length = stringHeader.get(index);

        if (stringBytes[index].hasArray()) {
            return new String(stringBytes[index].array(), 0, utf8length, "UTF-8");
        }
        byte[] bytes = new byte[utf8length];
        stringBytes[index].position(0);
        stringBytes[index].limit(utf8length +1 ); // account for extra zero
        stringBytes[index].get(bytes);

        return new String(bytes, 0, utf8length, "UTF-8");
    }

    public boolean getBoolean(int index) {
        byte rawByte = booleanBytes.get(index);

        return rawByte == TRUE_BYTE;

    }

    public int getInteger(int index) {
        return intBytes.getInt(index * SIZEOF_INT);
    }

    /**
     * Get all buffers, possibly including the buffers containing the strings.
     * 
     * @return all buffers.
     * 
     * @param includeStringBuffers
     *            if true, the buffers for holding the values of strings will be included.
     */
    public ByteBuffer[] getByteBuffers(boolean includeStringBuffers) {
        if (includeStringBuffers) {
            return byteBuffers;
        } else {
            return allButStringByteBuffers;
        }
    }

    public ByteBuffer[] getStringByteBuffers() {
        return stringBytes;
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
            throw new IOException("Amuse message in inconsistent state, strign count greater than number of string buffers");
        }

        for (int i = 0; i < getStringCount(); i++) {
            int utf8Length = stringHeader.get(i) +1; // account for extra zero

            stringBytes[i].clear().limit(utf8Length );
        }

        // set the limit of the rest of the string bytes to 0
        for (int i = getStringCount(); i < stringBytes.length; i++) {
            stringBytes[i].limit(0);
        }
    }

    public void writeTo(SocketChannel channel) throws IOException {
        //System.err.prinln("writing to socket channel: " + this.toContentString());
        //System.err.prinln("writing to socket channel: " + this);

        headerBytes.clear();
        setPrimitiveLimitsFromHeader();
        setStringLimitsFromHeader();

        // write to channel
        // channel.write(byteBuffers);

        // write all bufferd to channel
        boolean done = false;
        
        while(!done) {
            channel.write(byteBuffers);
            
            done = true;
            for (ByteBuffer buffer : byteBuffers) {
                if (buffer.hasRemaining()) {
                   done = false;
                }
            }
        }

        // alternative, debugging version of writing buffers
        // for (ByteBuffer buffer : byteBuffers) {
        // //System.err.println("writing " + buffer + " of length "
        // + buffer.remaining());
        // channel.write(buffer);
        //
        // if (buffer.hasRemaining()) {
        // System.err.println("Error! not all bytes written "
        // + buffer.remaining());
        // }
        // }
    }

    public void writeTo(WriteMessage writeMessage) throws IOException {
        if (logger.isTraceEnabled()) {
            logger.trace("writing to write message: {} ", this);
            //logger.trace("writing to write message: {} ", this.toContentString());
        }

        headerBytes.clear();
        setPrimitiveLimitsFromHeader();
        setStringLimitsFromHeader();

        for (ByteBuffer buffer : byteBuffers) {
            if (buffer.hasRemaining()) {
                writeMessage.writeByteBuffer(buffer);
            }
        }
    }

    // make sure there is enough space for each primitive buffer
    // (including the string header)
    public boolean ensurePrimitiveCapacity() {
        boolean buffersUpdated = false;

        if (getIntCount() * SIZEOF_INT > intBytes.capacity()) {
            intBytes = ByteBuffer.allocateDirect(getIntCount() * SIZEOF_INT);
            intBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getLongCount() * SIZEOF_LONG > longBytes.capacity()) {
            longBytes = ByteBuffer.allocateDirect(getLongCount() * SIZEOF_LONG);
            longBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getFloatCount() * SIZEOF_FLOAT > floatBytes.capacity()) {
            floatBytes = ByteBuffer.allocateDirect(getFloatCount() * SIZEOF_FLOAT);
            floatBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getDoubleCount() * SIZEOF_DOUBLE > doubleBytes.capacity()) {
            doubleBytes = ByteBuffer.allocateDirect(getDoubleCount() * SIZEOF_DOUBLE);
            doubleBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getBooleanCount() * SIZEOF_BOOLEAN > booleanBytes.capacity()) {
            booleanBytes = ByteBuffer.allocateDirect(getBooleanCount() * SIZEOF_BOOLEAN);
            booleanBytes.order(getByteOrder());
            buffersUpdated = true;
        }

        if (getStringCount() * SIZEOF_INT > stringHeaderBytes.capacity()) {
            stringHeaderBytes = ByteBuffer.allocateDirect(getStringCount() * SIZEOF_INT);
            stringHeaderBytes.order(getByteOrder());
            stringHeader = stringHeaderBytes.asIntBuffer();
            buffersUpdated = true;
        }

        if (buffersUpdated) {
            allButStringByteBuffers = new ByteBuffer[] { headerBytes, intBytes, longBytes, floatBytes, doubleBytes, booleanBytes,
                    stringHeaderBytes };

            // update byte buffers array
            ByteBuffer[] newByteBuffers = new ByteBuffer[allButStringByteBuffers.length + stringBytes.length];
            for (int i = 0; i < allButStringByteBuffers.length; i++) {
                newByteBuffers[i] = allButStringByteBuffers[i];
            }
            for (int i = 0; i < stringBytes.length; i++) {
                newByteBuffers[allButStringByteBuffers.length + i] = stringBytes[i];
            }
            byteBuffers = newByteBuffers;

            //System.err.println("ensurePrimitiveCapacity() Updated buffers to " + Arrays.toString(byteBuffers));
        }

        return buffersUpdated;
    }

    public boolean ensureStringsCapacity() {
        // checking if the string header is big enough is checked above, so
        // we
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
            int stringLength = stringHeader.get(i) +1 ; // account for extra zero
            if (stringBytes[i] == null || stringLength > stringBytes[i].capacity()) {

                stringBytes[i] = ByteBuffer.allocateDirect(stringLength);
                buffersUpdated = true;
            }
        }

        if (buffersUpdated) {
            // update byte buffers array
            ByteBuffer[] newByteBuffers = new ByteBuffer[allButStringByteBuffers.length + stringBytes.length];
            for (int i = 0; i < allButStringByteBuffers.length; i++) {
                newByteBuffers[i] = allButStringByteBuffers[i];
            }
            for (int i = 0; i < stringBytes.length; i++) {
                newByteBuffers[allButStringByteBuffers.length + i] = stringBytes[i];
            }
            byteBuffers = newByteBuffers;

            //System.err.println("ensureStringsCapacity() Updated buffers to " + Arrays.toString(byteBuffers));
        }

        return buffersUpdated;
    }

    public boolean ensureStringsCapacity(int index) {
        // checking if the string header is big enough is checked above, so
        // we
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

        if (buffersUpdated) {
            // update byte buffers array
            ByteBuffer[] newByteBuffers = new ByteBuffer[allButStringByteBuffers.length + stringBytes.length];
            for (int i = 0; i < allButStringByteBuffers.length; i++) {
                newByteBuffers[i] = allButStringByteBuffers[i];
            }
            for (int i = 0; i < stringBytes.length; i++) {
                newByteBuffers[allButStringByteBuffers.length + i] = stringBytes[i];
            }
            byteBuffers = newByteBuffers;

            
            //System.err.println("ensureStringsCapacity() Updated buffers to " + Arrays.toString(byteBuffers));
        }

        {
          int stringLength = stringHeader.get(index) +1; // account for extra zero
          if (stringBytes[index] == null || stringLength > stringBytes[index].capacity()) {

              stringBytes[index] = ByteBuffer.allocateDirect(stringLength);
              byteBuffers[allButStringByteBuffers.length + index] = stringBytes[index];
          }
        }


        return buffersUpdated;
    }

    public boolean readFrom(SocketChannel channel) throws IOException {
        boolean updatedBuffers = false;

        //System.err.println("receiving header from channel");

        headerBytes.clear();

        readAll(channel, headerBytes);

        // set buffers to byte order specified in buffer
        setByteOrder(getByteOrder());

        //System.err.println("reading content for " + this);

        if (ensurePrimitiveCapacity()) {
            updatedBuffers = true;
        }

        // then, set limits for primitive buffers, and receive those

        setPrimitiveLimitsFromHeader();

        //System.err.println("receiving primitives from channel");

        headerBytes.position(headerBytes.limit());
        // we also request to read the header, but its position is already
        // equal to its limit, so no bytes are read into it.
        readAll(channel, allButStringByteBuffers);

        // make sure there is enough space for the strings
        if (ensureStringsCapacity()) {
            updatedBuffers = true;
        }

        // set the limits
        setStringLimitsFromHeader();

        //System.err.println("receiving strings from channel");

        // and receive!
        readAll(channel, stringBytes);

        //System.err.println("done receiving message from channel: " + this);

        return updatedBuffers;
    }

    public boolean readFrom(ReadMessage readMessage) throws IOException {
        boolean updatedBuffers = false;

        logger.trace("reading header from message");

        headerBytes.clear();

        readMessage.readByteBuffer(headerBytes);

        // set buffers to byte order specified in buffer
        setByteOrder(getByteOrder());

        if (ensurePrimitiveCapacity()) {
            updatedBuffers = true;
        }

        // then, set limits for primitive buffers, and receive those

        setPrimitiveLimitsFromHeader();

        logger.trace("reading primitives from message");

        headerBytes.position(headerBytes.capacity());
        for (ByteBuffer buffer : allButStringByteBuffers) {
            if (buffer.hasRemaining()) {
                readMessage.readByteBuffer(buffer);
            }
        }

        // make sure there is enough space for the strings
        if (ensureStringsCapacity()) {
            updatedBuffers = true;
        }

        // set the limits
        setStringLimitsFromHeader();

        logger.trace("reading strings from message");

        // and receive!
        for (ByteBuffer buffer : stringBytes) {
            if (buffer.hasRemaining()) {
                readMessage.readByteBuffer(buffer);
            }
        }

        if (logger.isTraceEnabled()) {
            logger.trace("done receiving message from ReadMessage: " + this);
        }

        return updatedBuffers;
    }

    public String toContentString() throws IOException {
        String message = "AmuseMessage <id:" + getCallID() + " function ID:" + getFunctionID() + " count:" + getCallCount();

        if (isErrorState()) {
            message = message + " ERROR";
        }

        if (getByteOrder() == ByteOrder.BIG_ENDIAN) {
            message = message + " order: B";
        } else {
            message = message + " order: l";
        }

        if (getIntCount() != 0) {
            intBytes.clear();
            message = message + " ints: [";
            for (int i = 0; i < getIntCount(); i++) {
                message = message + ", " + intBytes.getInt(i * SIZEOF_INT);
            }
            message = message + "] ";
        }

        if (getLongCount() != 0) {
            longBytes.clear();
            message = message + " longs: [";
            for (int i = 0; i < getLongCount(); i++) {
                message = message + ", " + longBytes.getLong(i * SIZEOF_LONG);
            }
            message = message + "] ";
        }

        if (getFloatCount() != 0) {
            floatBytes.clear();
            message = message + " floats: [";
            for (int i = 0; i < getFloatCount(); i++) {
                message = message + ", " + floatBytes.getFloat(i * SIZEOF_FLOAT);
            }
            message = message + "] ";
        }

        if (getDoubleCount() != 0) {
            doubleBytes.clear();
            message = message + " double: [";
            for (int i = 0; i < getDoubleCount(); i++) {
                message = message + ", " + doubleBytes.getDouble(i * SIZEOF_DOUBLE);
            }
            message = message + "] ";
        }

        if (getBooleanCount() != 0) {
            message = message + " boolean: [";
            for (int i = 0; i < getBooleanCount(); i++) {
                message = message + ", " + getBoolean(i);
            }
            message = message + "] ";
        }

        if (getStringCount() != 0) {
            message = message + " string: [";
            for (int i = 0; i < getStringCount(); i++) {
                message = message + ", " + getString(i);
            }
            message = message + "] ";
        }

        message = message + ">";

        // return "Call <id:" + getCallID() + " function ID:" +
        // getFunctionID()
        // + " count:" + getCount() + " ints:" + getIntCount()
        // + " longs: " + getLongCount() + " floats:" + getFloatCount()
        // + " doubles:" + getDoubleCount() + " booleans:"
        // + getBooleanCount() + " strings:" + getStringCount()
        // + " byte order:" + getByteOrder() + " error:"
        // + isErrorState() + ">";

        return message;
    }

    public String toString() {
        String message = "AmuseMessage <id:" + getCallID() + " function ID:" + getFunctionID() + " count:" + getCallCount();

        if (isErrorState()) {
            message = message + " ERROR";
        }

        if (getByteOrder() == ByteOrder.BIG_ENDIAN) {
            message = message + " order: B";
        } else {
            message = message + " order: l";
        }

        if (getIntCount() != 0) {
            message = message + " ints:" + getIntCount();
        }

        if (getLongCount() != 0) {
            message = message + " longs:" + getLongCount();
        }

        if (getFloatCount() != 0) {
            message = message + " floats:" + getFloatCount();
        }

        if (getDoubleCount() != 0) {
            message = message + " doubles:" + getDoubleCount();
        }

        if (getBooleanCount() != 0) {
            message = message + " booleans:" + getBooleanCount();
        }

        if (getStringCount() != 0) {
            message = message + " strings:" + getStringCount();
        }

        message = message + ">";

        // return "Call <id:" + getCallID() + " function ID:" +
        // getFunctionID()
        // + " count:" + getCount() + " ints:" + getIntCount()
        // + " longs: " + getLongCount() + " floats:" + getFloatCount()
        // + " doubles:" + getDoubleCount() + " booleans:"
        // + getBooleanCount() + " strings:" + getStringCount()
        // + " byte order:" + getByteOrder() + " error:"
        // + isErrorState() + ">";

        return message;
    }

    public void setIntCount(int ints) {
        header.put(HEADER_INT_COUNT_INDEX, ints);
    }

    public void setLongCount(int longs) {
        header.put(HEADER_LONG_COUNT_INDEX, longs);
    }

    public void setFloatCount(int floats) {
        header.put(HEADER_FLOAT_COUNT_INDEX, floats);
    }

    public void setDoubleCount(int doubles) {
        header.put(HEADER_DOUBLE_COUNT_INDEX, doubles);
    }

    public void setBooleanCount(int booleans) {
        header.put(HEADER_BOOLEAN_COUNT_INDEX, booleans);
    }

    public void setStringCount(int strings) {
        header.put(HEADER_STRING_COUNT_INDEX, strings);
    }

    public int[] getIntSlice(int sliceIndex) {
        int[] result = new int[getCallCount()];

        intBytes.position(getCallCount() * sliceIndex * SIZEOF_INT);
        intBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_INT);

        intBytes.asIntBuffer().get(result);

        return result;
    }

    public long[] getLongSlice(int sliceIndex) {
        long[] result = new long[getCallCount()];

        longBytes.position(getCallCount() * sliceIndex * SIZEOF_LONG);
        longBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_LONG);

        longBytes.asLongBuffer().get(result);

        return result;
    }

    public float[] getFloatSlice(int sliceIndex) {
        float[] result = new float[getCallCount()];

        floatBytes.position(getCallCount() * sliceIndex * SIZEOF_FLOAT);
        floatBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_FLOAT);

        floatBytes.asFloatBuffer().get(result);

        return result;
    }

    public double[] getDoubleSlice(int sliceIndex) {
        double[] result = new double[getCallCount()];

        doubleBytes.position(getCallCount() * sliceIndex * SIZEOF_DOUBLE);
        doubleBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_DOUBLE);

        doubleBytes.asDoubleBuffer().get(result);

        return result;
    }

    public boolean[] getBooleanSlice(int sliceIndex) throws IOException {
        int callCount = getCallCount();

        boolean[] result = new boolean[callCount];

        int offset = sliceIndex * callCount;

        for (int i = 0; i < callCount; i++) {
            result[i] = getBoolean(offset + i);
        }

        return result;
    }

    public String[] getStringSlice(int sliceIndex) throws IOException {
        int callCount = getCallCount();

        String[] result = new String[callCount];

        int offset = sliceIndex * callCount;

        for (int i = 0; i < callCount; i++) {
            result[i] = getString(offset + i);
        }

        return result;
    }

    // sets all elements of a slice
    public void setIntSlice(int sliceIndex, int[] data) {
        intBytes.position(getCallCount() * sliceIndex * SIZEOF_INT);
        intBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_INT);

        intBytes.asIntBuffer().put(data);
    }

    // sets all elements of a slice to a single value
    public void setIntSlice(int sliceIndex, int value) {
        intBytes.position(getCallCount() * sliceIndex * SIZEOF_INT);
        intBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_INT);

        IntBuffer buffer = intBytes.asIntBuffer();

        while (buffer.hasRemaining()) {
            buffer.put(value);
        }
    }

    // sets all elements of a slice
    public void setLongSlice(int sliceIndex, long[] data) {
        longBytes.position(getCallCount() * sliceIndex * SIZEOF_LONG);
        longBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_LONG);

        longBytes.asLongBuffer().put(data);
    }

    // sets all elements of a slice to a single value
    public void setLongSlice(int sliceIndex, long value) {
        longBytes.position(getCallCount() * sliceIndex * SIZEOF_LONG);
        longBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_LONG);

        LongBuffer buffer = longBytes.asLongBuffer();

        while (buffer.hasRemaining()) {
            buffer.put(value);
        }
    }

    // sets all elements of a slice
    public void setFloatSlice(int sliceIndex, float[] data) {
        floatBytes.position(getCallCount() * sliceIndex * SIZEOF_FLOAT);
        floatBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_FLOAT);

        floatBytes.asFloatBuffer().put(data);
    }

    // sets all elements of a slice to a single value
    public void setFloatSlice(int sliceIndex, float value) {
        floatBytes.position(getCallCount() * sliceIndex * SIZEOF_FLOAT);
        floatBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_FLOAT);

        FloatBuffer buffer = floatBytes.asFloatBuffer();

        while (buffer.hasRemaining()) {
            buffer.put(value);
        }
    }

    // sets all elements of a slice
    public void setDoubleSlice(int sliceIndex, double[] data) {
        doubleBytes.position(getCallCount() * sliceIndex * SIZEOF_DOUBLE);
        doubleBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_DOUBLE);

        doubleBytes.asDoubleBuffer().put(data);
    }

    // sets all elements of a slice to a single value
    public void setDoubleSlice(int sliceIndex, double value) {
        doubleBytes.position(getCallCount() * sliceIndex * SIZEOF_DOUBLE);
        doubleBytes.limit(getCallCount() * (sliceIndex + 1) * SIZEOF_DOUBLE);

        DoubleBuffer buffer = doubleBytes.asDoubleBuffer();

        while (buffer.hasRemaining()) {
            buffer.put(value);
        }
    }

    // sets all elements of a slice
    public void setBooleanSlice(int sliceIndex, boolean[] data) {
        int callCount = getCallCount();
        for (int i = 0; i < callCount; i++) {
            setBoolean((callCount * sliceIndex) + i, data[i]);
        }
    }

    // sets all elements of a slice to a single value
    public void setBooleanSlice(int sliceIndex, boolean value) {
        int callCount = getCallCount();
        for (int i = 0; i < callCount; i++) {
            setBoolean((callCount * sliceIndex) + i, value);
        }
    }

    // sets all elements of a slice
    public void setStringSlice(int sliceIndex, String[] data) {
        int callCount = getCallCount();
        for (int i = 0; i < callCount; i++) {
            setString((callCount * sliceIndex) + i, data[i]);
        }
    }

    // sets all elements of a slice to a single value
    public void setStringSlice(int sliceIndex, String value) {
        int callCount = getCallCount();
        for (int i = 0; i < callCount; i++) {
            setString((callCount * sliceIndex) + i, value);
        }
    }

}
