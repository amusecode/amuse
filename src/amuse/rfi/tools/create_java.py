from amuse.support.core import late
from amuse.support import exceptions
from amuse import config
from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.core import LegacyFunctionSpecification

import sys
import os
import inspect

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('Int', 'Int', '', 'int', ''),
    'int64' : DTypeSpec('Long', 'Long',
                    '', 'long', ''),
    'float32' : DTypeSpec('Float', 'Float',
                    '', 'float', ''),
    'float64' : DTypeSpec('Double', 'Double',
                    '', 'double', ''),
    'bool' : DTypeSpec('Boolean', 'Boolean',
                    '', 'boolean', ''),
    'string' : DTypeSpec('String', 'String',
                    '', 'String', ''),
})

IMPORTS_CODE_STRING = """
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.InetSocketAddress;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.nio.FloatBuffer;
import java.nio.DoubleBuffer;
import java.nio.channels.SocketChannel;
import java.util.Arrays;
"""

AMUSE_MESSAGE_CLASS_CODE_STRING = """
 private static class AmuseMessage {
        public static final int HEADER_SIZE = 10; // integers

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

        private static final long serialVersionUID = 1L;

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

        AmuseMessage(int callID, int functionID, int count) {
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
        AmuseMessage(int callID, int functionID, int count, String error) {
            this();

            setCallID(callID);
            setFunctionID(functionID);
            setCallCount(count);
            setError(error);
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
         * Change the byte order of this message. Also swaps the content of all
         * the buffers, if requested.
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

            // set error state
            headerBytes.put(HEADER_ERROR_FLAG, TRUE_BYTE);

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
                bytes = value.getBytes("UTF-8");

                // set length of string in header
                stringHeader.put(position, bytes.length);

                // make sure there is space for the string
                ensureStringsCapacity();

                stringBytes[position].clear();
                stringBytes[position].put(bytes);

            } catch (UnsupportedEncodingException e) {
                System.err.println("ERROR! UTF-8 not supported by the JVM!");
            }
        }
        
        public void setString(int index, String value) {
            // encode string to UTF-8
            byte[] bytes;

            try {
                bytes = value.getBytes("UTF-8");

                // set length of string in header
                stringHeader.put(index, bytes.length);

                // make sure there is space for the string
                ensureStringsCapacity();

                stringBytes[index].clear();
                stringBytes[index].put(bytes);

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
            stringBytes[index].limit(utf8length);
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
         * Get all buffers, possibly including the buffers containing the
         * strings.
         * 
         * @return all buffers.
         * 
         * @param includeStringBuffers
         *            if true, the buffers for holding the values of strings
         *            will be included.
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
            //System.err.prinln("writing to socket channel: " + this.toContentString());
            //System.err.prinln("writing to socket channel: " + this);
            

            headerBytes.clear();
            setPrimitiveLimitsFromHeader();
            setStringLimitsFromHeader();

            // write to channel
            channel.write(byteBuffers);

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
                allButStringByteBuffers = new ByteBuffer[] { headerBytes, intBytes, longBytes, floatBytes, doubleBytes,
                        booleanBytes, stringHeaderBytes };

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
                int stringLength = stringHeader.get(i);
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

        boolean readFrom(SocketChannel channel) throws IOException {
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

        public String toContentString() throws IOException {
            String message = "AmuseMessage <id:" + getCallID() + " function ID:" + getFunctionID() + " count:"
                    + getCallCount();

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
            String message = "AmuseMessage <id:" + getCallID() + " function ID:" + getFunctionID() + " count:"
                    + getCallCount();

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
            
            for(int i = 0; i < callCount; i++) {
                result[i] = getBoolean(offset + i);
            }
            
            return result;
        }
        
                
        public String[] getStringSlice(int sliceIndex) throws IOException {
            int callCount = getCallCount();
         
            String[] result = new String[callCount];

            int offset = sliceIndex * callCount;
            
            for(int i = 0; i < callCount; i++) {
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
            
            while(buffer.hasRemaining()) {
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
            
            while(buffer.hasRemaining()) {
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
            
            while(buffer.hasRemaining()) {
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
            
            while(buffer.hasRemaining()) {
                buffer.put(value);
            }
        }

       // sets all elements of a slice
       public void setBooleanSlice(int sliceIndex, boolean[] data) {
           int callCount = getCallCount();
           for(int i = 0; i < callCount; i++) {
               setBoolean((callCount * sliceIndex) + i, data[i]);
            }
        }
        
        // sets all elements of a slice to a single value
        public void setBooleanSlice(int sliceIndex, boolean value) {
           int callCount = getCallCount();
           for(int i = 0; i < callCount; i++) {
               setBoolean((callCount * sliceIndex) + i, value);
            }
        }
        
        // sets all elements of a slice
        public void setStringSlice(int sliceIndex, String[] data) {
           int callCount = getCallCount();
           for(int i = 0; i < callCount; i++) {
               setString((callCount * sliceIndex) + i, data[i]);
            }
        }
        
        // sets all elements of a slice to a single value
        public void setStringSlice(int sliceIndex, String value) {
           int callCount = getCallCount();
           for(int i = 0; i < callCount; i++) {
               setString((callCount * sliceIndex) + i, value);
            }
        }
        
    }
"""



FOOTER_CODE_STRING = """
  private final AmuseMessage request;
  private final AmuseMessage reply;
  private final CodeInterface code;
    
  Worker() {
      this.request = new AmuseMessage();
      this.reply = new AmuseMessage();
      
      code = new Code();
  }

  private void runSockets(int port) {
        try {
            SocketChannel channel = SocketChannel.open(new InetSocketAddress(port));
            channel.socket().setTcpNoDelay(true);

            boolean keepRunning = true;
            while (keepRunning) {
                request.clear();
                request.readFrom(channel);

                //System.err.println("got message " + request.toString());

                reply.clear();

                reply.setCallID(request.getCallID());
                reply.setFunctionID(request.getFunctionID());
                reply.setCallCount(request.getCallCount());

                keepRunning = handleCall();

                //System.err.println("sending reply message " + reply.toString());
                //System.err.println("sending reply message " + reply.toContentString());

                reply.writeTo(channel);
                
                //System.err.println("call handled");
            }
        } catch (IOException e) {
            System.err.println("Error running worker: " + e.getMessage());
        }
    }

    public static void main(String[] arguments) throws IOException {
        //System.err.println("Java worker");
        //for (String argument : arguments) {
        //    System.err.println("argument: " + argument);
        //}

        if (arguments.length == 0) {
            System.err.println("No arguments to java worker. expected a socket port number");
            System.exit(1);
        }

        int port = Integer.parseInt(arguments[0]);

        new Worker().runSockets(port);

    }    



"""

class MakeJavaCodeString(GenerateASourcecodeString):
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
       
         

class GenerateAJavaStringOfAFunctionSpecification(MakeJavaCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
   
        
    def start(self):
        #must and can handle array is the same thing in Java codes...
        if self.specification.can_handle_array:
            self.specification.must_handle_array = True

        self.specification.prepare_output_parameters()
        self.output_casestmt_start()
        self.out.indent()
        self.out.lf() + "{"
        self.out.indent()
        
        self.output_lines_with_number_of_outputs()
        
        
        if self.specification.name.startswith("internal__"):
            self.out.lf() + "//" +  self.specification.name + " ignored"
        else:
            self.output_declare_variables()
            self.output_function_start()
            self.output_function_parameters()
            self.output_function_end()
            self.output_copy_output_variables()

        self.out.dedent()
        self.out.lf() + "}"
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
        
    def output_casestmt_start(self):
        self.out + 'case ' + self.specification.id + ':'
        
              
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for parameter in self.specification.output_parameters:
            count = dtype_to_count.get(parameter.datatype, 0)
            dtype_to_count[parameter.datatype] = count + 1
                
        if not self.specification.result_type is None:
            count = dtype_to_count.get(self.specification.result_type, 0)
            dtype_to_count[self.specification.result_type] = count + 1
            
        for dtype in dtype_to_count:       
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_count[dtype]
            self.out.lf() + 'reply.set' + spec.input_var_name + 'Count(' + count + ' * count);'
            pass
        
        self.out.lf() + 'reply.ensurePrimitiveCapacity();'
            
    
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ', '
                
            if parameter.direction == LegacyFunctionSpecification.IN:
                if self.specification.must_handle_array:
                    self.out + parameter.name
                else:
                    self.out + parameter.name + '[0]'
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                    self.out + parameter.name
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                    self.out + parameter.name
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out + 'count'
    
        self.out.dedent()

    def output_declare_variables(self):
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out.lf() + spec.type + ' result;'
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.IN or parameter.direction == LegacyFunctionSpecification.INOUT :
                self.out.lf() + spec.type + '[] ' + parameter.name + ' = request.get' + spec.input_var_name + 'Slice(' + parameter.input_index + ');'
            if parameter.direction == LegacyFunctionSpecification.OUT:
                self.out.lf() + spec.type + '[] ' + parameter.name + ' = new ' +  spec.type + '[count];'

    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            self.out + 'result = '
            
        self.out + 'code.' + self.specification.name + '('

    def output_function_end(self):
        self.out + ')' + ';'
                
    def output_copy_output_variables(self):
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out.lf() + 'reply.set' + spec.output_var_name + 'Slice(0, result);'

        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.OUT or parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() + 'reply.set' + spec.output_var_name + 'Slice(' + parameter.output_index + ', ' + parameter.name + ');'
  
    def output_casestmt_end(self):
        self.out.n() + 'break;'

class GenerateAJavaFunctionDeclarationStringFromAFunctionSpecification(MakeJavaCodeString):
   
        
    def start(self):
        #must and can handle array is the same thing in Java codes...
        if self.specification.can_handle_array:
            self.specification.must_handle_array = True
             
        self.output_function_parameter_types()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self._result = self.out.string
        
    def output_function_parameter_types(self):        
        for parameter in self.specification.parameters:
            if (parameter.direction == LegacyFunctionSpecification.IN):
                self.out.lf() + '// parameter "' + parameter.name + '" is an input parameter'
            elif (parameter.direction == LegacyFunctionSpecification.OUT):
                self.out.lf() + '// parameter "' + parameter.name + '" is an output parameter'
            elif (parameter.direction == LegacyFunctionSpecification.INOUT):
                self.out.lf() + '// parameter "' + parameter.name + '" is an inout parameter'
            elif (parameter.direction == LegacyFunctionSpecification.LENGTH):
                self.out.lf() + '// parameter "' + parameter.name + '" is a length parameter'
            
    def output_function_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ', '
                
            self.out + spec.type
            if ((self.specification.must_handle_array and parameter.is_input()) or parameter.is_output()):
                self.out + '[]'
                
            self.out + ' '
            self.out + parameter.name
                
            
    def output_function_end(self):
        self.out + ')' + ';'
        
    def output_function_start(self):
        self.out.n()
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.type
            self.out + ' '
        else:
            self.out + 'void' + ' '
        self.out + self.specification.name + '('
        
class GenerateAJavaSourcecodeStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def specification_class(self):
        raise exceptions.AmuseException("No specification_class set, please set the specification_class first")
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def output_sourcecode_for_function(self):
        return GenerateAJavaStringOfAFunctionSpecification()
    
    def start(self):
        
        self.out.lf()

        self.out + IMPORTS_CODE_STRING
        
        self.output_imports()
        
        self.out.lf() + 'class Worker {'
        self.out.indent().lf()
        self.out + AMUSE_MESSAGE_CLASS_CODE_STRING

        
        self.output_handle_call()
        
        self.out.lf() + FOOTER_CODE_STRING

        self.out.dedent().lf()
        
        self.out.lf() + "}"
        
        self._result = self.out.string
        
    def output_imports(self):
        if hasattr(self.specification_class, 'imports'):
            for clazz in self.specification_class.imports:
                self.out.n() + 'import ' + clazz + ';'
        self.out.lf()
    
    def output_code_constants(self):
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            
            maxin = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
            self.out + 'static int MAX_' + dtype_spec.input_var_name.upper() + ' = ' + maxin + ";"
            self.out.lf()
            
            maxout = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            self.out + 'static int MAX_' + dtype_spec.output_var_name.upper() + ' = ' + maxout + ";"
            self.out.lf()
            
    def output_handle_call(self):
        
        self.out.lf().lf() + 'private boolean handleCall() throws IOException {'
        self.out.indent()
        
        self.out.lf() + 'int count = request.getCallCount();'
        
        self.out.lf().lf() + 'switch (request.getFunctionID()) {'
        self.out.indent()
        self.out.lf() + 'case 0:'
        self.out.indent().lf() + 'return false;'
        self.out.dedent()
        
        self.output_sourcecode_for_functions()
        
        self.out.lf() + 'default:'
        self.out.indent()
        self.out.lf() + 'System.err.println("unknown function id " + request.getFunctionID());'
        self.out.lf() + 'reply.setError("unknown function id " + request.getFunctionID());'
        self.out.dedent()
        
        self.out.dedent().lf() + '}'
        self.out.dedent()
        self.out.indent().lf() + 'return true;'
        self.out.dedent().lf() + '}'

class GenerateAJavaInterfaceStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def ignore_functions_from_specification_classes(self):
        return []
        
    @late
    def underscore_functions_from_specification_classes(self):
        return []
        
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
    
        
    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        for cls in self.ignore_functions_from_specification_classes:
            if hasattr(cls, x.specification.name):
                return False
        
        return True
        
    def output_sourcecode_for_function(self):
        return GenerateAJavaFunctionDeclarationStringFromAFunctionSpecification()
        
    def start(self):  
        self.out + 'public interface CodeInterface {'
        self.out.indent().lf()
            
        self.output_sourcecode_for_functions()
        
        self.out.dedent().lf() + '}'
        
        self.out.lf()
        
        self._result = self.out.string
     
class GenerateAJavaWorkerScript(GenerateASourcecodeString):

    @late
    def code_dir(self):
        return os.getcwd()

    @late
    def java(self):
        return config.java.java

    @late
    def template_dir(self):
        return os.path.dirname(__file__)
        
    @late
    def template_string(self):
        path = self.template_dir
        path = os.path.join(path, 'java_code_script.template')
            
        with open(path, "r") as f:
            template_string = f.read()
        
        return template_string
    
    def script_string(self):
        return self.template_string.format(
            executable = sys.executable,
            code_dir = self.code_dir,
            java = self.java,
            classpath = self.specification_class.classpath
            )


    #def code_directory(self):
        #interface_module = inspect.getmodule(self.specification_class).__name__
        #return os.path.dirname(inspect.getfile(self.specification_class))
    
    def start(self):
        self.out + self.script_string()
        
        self._result = self.out.string

