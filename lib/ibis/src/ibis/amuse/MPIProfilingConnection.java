package ibis.amuse;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MPIProfilingConnection extends Thread {

	private static final Logger logger = LoggerFactory
			.getLogger(MPIProfilingConnection.class);

	public static final int SIZEOF_INT = 4;
	public static final int SIZEOF_LONG = 8;
	public static final int SIZEOF_FLOAT = 4;
	public static final int SIZEOF_DOUBLE = 8;
	public static final int SIZEOF_BOOLEAN = 1;
	
	public static final int MAGIC_NUMBER = 449682;

	private final SocketChannel channel;
	private final MPIProfilingCollector collector;
	private final int poolSize;

	public MPIProfilingConnection(SocketChannel channel,
			MPIProfilingCollector collector, int poolSize) {
		this.channel = channel;
		this.collector = collector;
		this.poolSize = poolSize;

		this.setDaemon(true);
		this.setName(this.getClass().getName());
		this.start();
	}

	void readAll(ByteBuffer buffer) throws IOException {
		while (buffer.hasRemaining()) {
			long read = channel.read(buffer);

			if (read == -1) {
				throw new IOException("Connection closed on reading data");
			}
		}
	}

	public void run() {
		try {

			ByteBuffer buffer = ByteBuffer.allocateDirect((SIZEOF_INT * 2) + (poolSize * SIZEOF_LONG));
			buffer.order(ByteOrder.nativeOrder());

			buffer.clear().limit(SIZEOF_INT);

			readAll(buffer);
			
			if (buffer.getInt(0) != MAGIC_NUMBER) {
				throw new Exception("First int of message " + buffer.getInt(0) + " not equal to magic number " + MAGIC_NUMBER);
			}

			while (channel.isConnected()) {
				buffer.clear();
				readAll(buffer);
				
				buffer.clear();
				
				int rank = buffer.getInt();
				int size = buffer.getInt();
				
				long[] sentPerRank = new long[size];
				
				buffer.asLongBuffer().get(sentPerRank);
				
				collector.addSentBytes(rank, size, sentPerRank);
			}

		} catch (Exception e) {
			logger.error("Error on handling mpi profiling connection", e);

		}
	}

}
