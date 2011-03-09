package ibis.amuse;

import java.io.IOException;
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.nio.channels.SocketChannel;
import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Call implements Serializable {

    private static final long serialVersionUID = 1L;

    private static final Logger logger = LoggerFactory.getLogger(Call.class);

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


    private final int id;
    
    private final int functionID;
    
    // number of function calls in one "call message"
    private final int count;

    private final int[] ints;

    private final long[] longs;

    private final float[] floats;

    private final double[] doubles;

    private final boolean[] booleans;

    private final String[] strings;

    /**
     * Read a call from a socket
     * 
     * @throws IOException
     *             in case the call cannot be read
     */
    Call(SocketChannel channel) throws IOException {
        ByteBuffer headerBytes = ByteBuffer.allocate(MessageFormat.HEADER_SIZE * MessageFormat.SIZEOF_INT);
        headerBytes.order(ByteOrder.nativeOrder());

        logger.debug("reading header");
        readAll(channel, headerBytes);

        headerBytes.clear();
        IntBuffer header = headerBytes.asIntBuffer();

        id = header.get();
        functionID = header.get();
        count = header.get();

        // create arrays
        headerBytes.clear();
        ints = new int[header.get()];
        longs = new long[header.get()];
        floats = new float[header.get()];
        doubles = new double[header.get()];
        booleans = new boolean[header.get()];
        strings = new String[header.get()];

        ByteBuffer intBytes = ByteBuffer.allocate(ints.length * MessageFormat.SIZEOF_INT);
        ByteBuffer longBytes = ByteBuffer.allocate(longs.length * MessageFormat.SIZEOF_LONG);
        ByteBuffer floatBytes = ByteBuffer.allocate(floats.length
                * MessageFormat.SIZEOF_FLOAT);
        ByteBuffer doubleBytes = ByteBuffer.allocate(doubles.length
                * MessageFormat.SIZEOF_DOUBLE);
        // booleans as send as bytes
        ByteBuffer booleanBytes = ByteBuffer.allocate(booleans.length
                * MessageFormat.SIZEOF_BOOLEAN);
        // strings are preceded by a header with integers denoting
        // the length of each string
        ByteBuffer stringHeaderBytes = ByteBuffer.allocate(strings.length
                * MessageFormat.SIZEOF_INT);

        logger.trace("reading data for call " + this);

        // receive all data
        readAll(channel, intBytes, longBytes, floatBytes, doubleBytes,
                booleanBytes, stringHeaderBytes);

        // set order of all buffers to native order
        intBytes.order(ByteOrder.nativeOrder());
        longBytes.order(ByteOrder.nativeOrder());
        floatBytes.order(ByteOrder.nativeOrder());
        doubleBytes.order(ByteOrder.nativeOrder());
        booleanBytes.order(ByteOrder.nativeOrder());
        stringHeaderBytes.order(ByteOrder.nativeOrder());

        // reset position for reading
        intBytes.clear();
        longBytes.clear();
        floatBytes.clear();
        doubleBytes.clear();
        booleanBytes.clear();
        stringHeaderBytes.clear();

        // read data into arrays
        intBytes.asIntBuffer().get(ints);
        longBytes.asLongBuffer().get(longs);
        floatBytes.asFloatBuffer().get(floats);
        doubleBytes.asDoubleBuffer().get(doubles);

        for (int i = 0; i < booleans.length; i++) {
            byte value = booleanBytes.get(i);
            booleans[i] = (value == MessageFormat.TRUE);
        }

        // retreive lengths (in bytes!) of strings from header, create buffers
        // with suitable lengths
        IntBuffer stringLengths = stringHeaderBytes.asIntBuffer();
        ByteBuffer[] stringByteBuffers = new ByteBuffer[strings.length];
        for (int i = 0; i < strings.length; i++) {
            stringByteBuffers[i] = ByteBuffer.allocate(stringLengths.get(i));
            stringByteBuffers[i].order(ByteOrder.nativeOrder());
        }

        // receive all strings (as bytes)
        readAll(channel, stringByteBuffers);

        for (int i = 0; i < strings.length; i++) {
            stringByteBuffers[i].clear();
            strings[i] = new String(stringByteBuffers[i].array(), "UTF-8");
        }
    }

    public int getId() {
        return id;
    }

    public int getFunctionID() {
        return functionID;
    }

    public int getCount() {
        return count;
    }

    public int[] getInts() {
        return ints;
    }

    public long[] getLongs() {
        return longs;
    }

    public float[] getFloats() {
        return floats;
    }

    public double[] getDoubles() {
        return doubles;
    }

    public boolean[] getBooleans() {
        return booleans;
    }

    public String[] getStrings() {
        return strings;
    }
    
    public String getString(int index) throws IOException {
        if (strings == null || (strings.length + 1) < index) {
            throw new IOException("cannot get string at index " + index + " in call" + this);
        }
        return strings[index];
    }

    public String valuesToString() {
        return Arrays.toString(ints) + Arrays.toString(longs)
                + Arrays.toString(floats) + Arrays.toString(doubles)
                + Arrays.toString(booleans) + Arrays.toString(strings);
    }

    public String toString() {
        return "Call: function ID:" + functionID + " count:" + count + " ints:"
                + ints.length + " longs: " + longs.length + " floats:"
                + floats.length + " doubles:" + doubles.length + " booleans:"
                + booleans.length + " strings:" + strings.length;
    }

    // private String readString() throws IOException {
    // int length = readInt();
    //
    // logger.debug("reading string of size " + length);
    //
    // ByteBuffer bytes = ByteBuffer.allocate(length);
    //
    // while (bytes.hasRemaining()) {
    // int read = in.read(bytes.array(), bytes.position(), bytes
    // .capacity());
    //
    // if (read == -1) {
    // throw new IOException("Connection closed on reading integer");
    // }
    // bytes.position(bytes.position() + read);
    // }
    //
    // return new String(bytes.array(), "UTF-8");
    // }
    //
    // private void writeInt(int integer) throws IOException {
    // ByteBuffer bytes = ByteBuffer.allocate(4);
    // bytes.putInt(integer);
    //
    // out.write(bytes.array());
    // }
    //
    // private void writeString(String string) throws IOException {
    //
    // byte[] bytes = string.getBytes("UTF-8");
    //
    // writeInt(bytes.length);
    // out.write(bytes);
    // }
    //
    // // big-endian bytes converted to an int
    // private int readInt() throws IOException {
    // ByteBuffer bytes = ByteBuffer.allocate(4);
    //
    // while (bytes.hasRemaining()) {
    // int read = in.read(bytes.array(), bytes.position(), bytes
    // .capacity());
    //
    // if (read == -1) {
    // throw new IOException("Connection closed on reading integer");
    // }
    // bytes.position(bytes.position() + read);
    // }
    //
    // return bytes.getInt(0);
    // }

}
