package ibis.amuse;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.net.Socket;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Call implements Serializable {

	private static final Logger logger = LoggerFactory.getLogger(Call.class);

	private static final int HEADER_SIZE = 8;

	private static final int SIZEOF_INT = 4;
	private static final int SIZEOF_DOUBLE = 4;
	private static final int SIZEOF_FLOAT = 8;
	// private static final int SIZEOF_STRING = 4;
	private static final int SIZEOF_BOOLEAN = 4;
	private static final int SIZEOF_LONG = 8;

	private static final int FUNCTION_ID_INDEX = 0;
	private static final int COUNT_INDEX = 1;
	private static final int INTEGER_INDEX = 2;
	private static final int LONG_INDEX = 3;
	private static final int FLOAT_INDEX = 4;
	private static final int DOUBLE_INDEX = 5;
	private static final int BOOLEAN_INDEX = 6;
	private static final int STRING_INDEX = 7;

	private static void readAll(ByteBuffer bytes, SocketChannel channel)
			throws IOException {
		while (bytes.hasRemaining()) {
			int read = channel.read(bytes);

			if (read == -1) {
				throw new IOException("Connection closed on reading data");
			}
		}
	}

	private static void readAll(ByteBuffer[] bytes, SocketChannel channel)
			throws IOException {
		
		while (true) {
			long read = channel.read(bytes);

			if (read == 0) {
				//all buffers now full
				return;
			}
			
			if (read == -1) {
				throw new IOException("Connection closed on reading data");
			}
		}
	}

	private final int functionID;

	// number of function calls in one "call message"
	private final int count;

	private final int[] integers;
	
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
		ByteBuffer header = ByteBuffer.allocate(HEADER_SIZE);
		header.order(ByteOrder.nativeOrder());

		readAll(header, channel);

		functionID = header.getInt(FUNCTION_ID_INDEX);
		count = header.getInt(COUNT_INDEX);

		//create arrays
		integers = new int[header.getInt(INTEGER_INDEX)];
		longs = new long[header.getInt(LONG_INDEX)];
		floats = new float[header.getInt(FLOAT_INDEX)];
		doubles = new double[header.getInt(DOUBLE_INDEX)];
		booleans = new boolean[header.getInt(BOOLEAN_INDEX)];
		strings = new String[header.getInt(STRING_INDEX)];
		
		ByteBuffer integerBuffer = ByteBuffer.allocate(integers.length * SIZEOF_INT);
		integerBuffer.order(ByteOrder.nativeOrder());

		// receive all data
		readAll(integerBuffer, channel);

		// read data into arrays
		integerBuffer.asIntBuffer().get(integers);

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
