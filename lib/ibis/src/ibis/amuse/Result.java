package ibis.amuse;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.nio.channels.SocketChannel;

public class Result {

    private final int callID;

    private final int functionID;

    // number of function calls in one "call message"
    private final int count;

    private final Throwable exception;

    private final int[] ints;

    private final long[] longs;

    private final float[] floats;

    private final double[] doubles;

    private final boolean[] booleans;

    private final String[] strings;

    /**
     * Empty (void) result. 
     * 
     * @param callID
     *            id of the call
     * @param functionID
     *            function id 
     * @param count
     *            number of calls done
     */
    Result(int callID, int functionID, int count) {
        this.callID = callID;
        this.functionID = functionID;
        this.count = count;
        this.exception = null;
        ints = new int[0];
        longs = new long[0];
        floats = new float[0];
        doubles = new double[0];
        booleans = new boolean[0];
        strings = new String[0];
    }
    
    /**
     * Result with an exception. Will add exception message as a string, and set
     * count to 0 to denote a failure.
     * 
     * @param callID
     *            id of the call that generated the exception
     * @param functionID
     *            function id tried
     * @param exception
     *            the exception that occured
     */
    Result(int callID, int functionID, Throwable exception) {
        this.callID = callID;
        this.functionID = functionID;
        this.count = 0;
        this.exception = exception;
        ints = new int[0];
        longs = new long[0];
        floats = new float[0];
        doubles = new double[0];
        booleans = new boolean[0];
        strings = new String[] { exception.getMessage() };
    }

    Result(int callID, int functionID, int count, int[] ints, long[] longs,
            float[] floats, double[] doubles, boolean[] booleans,
            String[] strings) {
        this.callID = callID;
        this.functionID = functionID;
        this.count = count;
        this.exception = null;
        this.ints = ints;
        this.longs = longs;
        this.floats = floats;
        this.doubles = doubles;
        this.booleans = booleans;
        this.strings = strings;
    }
    
    private static void writeAll(SocketChannel channel, ByteBuffer... buffers) throws IOException {
        channel.write(buffers);
    }

    public void writeTo(SocketChannel channel) throws IOException {
        ByteBuffer headerBytes = ByteBuffer.allocate(MessageFormat.HEADER_SIZE
                * MessageFormat.SIZEOF_INT);
        headerBytes.order(ByteOrder.nativeOrder());

        IntBuffer header = headerBytes.asIntBuffer();
        header.put(callID);
        header.put(functionID);
        header.put(count);
        header.put(ints.length);
        header.put(longs.length);
        header.put(floats.length);
        header.put(doubles.length);
        header.put(booleans.length);
        header.put(strings.length);

        ByteBuffer intBytes = ByteBuffer.allocate(ints.length
                * MessageFormat.SIZEOF_INT);
        intBytes.order(ByteOrder.nativeOrder());
        intBytes.asIntBuffer().put(ints);

        ByteBuffer longBytes = ByteBuffer.allocate(longs.length
                * MessageFormat.SIZEOF_LONG);
        longBytes.order(ByteOrder.nativeOrder());
        longBytes.asLongBuffer().put(longs);

        ByteBuffer floatBytes = ByteBuffer.allocate(floats.length
                * MessageFormat.SIZEOF_FLOAT);
        floatBytes.order(ByteOrder.nativeOrder());
        floatBytes.asFloatBuffer().put(floats);

        ByteBuffer doubleBytes = ByteBuffer.allocate(doubles.length
                * MessageFormat.SIZEOF_DOUBLE);
        doubleBytes.order(ByteOrder.nativeOrder());
        doubleBytes.asDoubleBuffer().put(doubles);

        ByteBuffer booleanBytes = ByteBuffer.allocate(doubles.length
                * MessageFormat.SIZEOF_BOOLEAN);
        doubleBytes.order(ByteOrder.nativeOrder());
        for (boolean value : booleans) {
            if (value) {
                booleanBytes.put(MessageFormat.TRUE);
            } else {
                booleanBytes.put(MessageFormat.FALSE);
            }
        }
        booleanBytes.clear();

        ByteBuffer stringHeaderBytes = ByteBuffer.allocate(strings.length
                * MessageFormat.SIZEOF_INT);
        stringHeaderBytes.order(ByteOrder.nativeOrder());
        IntBuffer stringHeader = stringHeaderBytes.asIntBuffer();
        
        ByteBuffer[] stringBytes = new ByteBuffer[strings.length];
        
        for(int i = 0; i < strings.length; i++) {
            byte[] bytes = strings[i].getBytes("UTF-8");
            stringHeader.put(bytes.length);
            stringBytes[i] = ByteBuffer.wrap(bytes);
        }
        
        writeAll(channel, headerBytes, intBytes, longBytes, floatBytes, doubleBytes, booleanBytes, stringHeaderBytes);
        writeAll(channel, stringBytes);
    }
}
