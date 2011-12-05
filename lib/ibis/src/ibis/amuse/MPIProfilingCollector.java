package ibis.amuse;

import ibis.ipl.IbisIdentifier;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.net.InetSocketAddress;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;
import java.util.HashMap;
import java.util.Map;

import javax.management.MBeanServer;
import javax.management.ObjectName;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MPIProfilingCollector extends Thread implements
		MPIProfilingCollectorMBean {

	private static final Logger logger = LoggerFactory
			.getLogger(MPIProfilingCollector.class);

        public static final int MPI_PROFILING_PORT = 52861;
	
	private final Map<IbisIdentifier, Long> sent;

	private final IbisIdentifier[] rankToIbis;

	// connection with codes
	private final ServerSocketChannel serverSocket;

	public MPIProfilingCollector(IbisIdentifier[] rankToIbis) throws Exception {
		this.rankToIbis = rankToIbis;

		sent = new HashMap<IbisIdentifier, Long>();

		MBeanServer mbs = ManagementFactory.getPlatformMBeanServer();
		ObjectName name = new ObjectName(
				"ibis.amuse:type=MPIProfilingCollector");
		mbs.registerMBean(this, name);

		serverSocket = ServerSocketChannel.open();
		serverSocket.socket().bind(new InetSocketAddress(MPI_PROFILING_PORT));

		this.setDaemon(true);
		this.setName(this.getClass().getName());
		this.start();
	}

	@Override
	public synchronized Map<IbisIdentifier, Long> getSentBytesPerIbis() {
		return new HashMap<IbisIdentifier, Long>(sent);
	}

	private synchronized void addSentBytes(IbisIdentifier ibis, long bytesSent) {
		Long currentValue = sent.get(ibis);
		if (currentValue == null) {
			currentValue = 0L;
		}

		sent.put(ibis, currentValue + bytesSent);
		
		if (logger.isDebugEnabled()) {
		    logger.debug("ibis " + ibis + " now has " + (currentValue + bytesSent) + " sent bytes");
		}
	}

	synchronized void addSentBytes(int rank, int size, long[] bytesSentPerRank) {
		for(int i = 0; i < bytesSentPerRank.length; i++) {
			addSentBytes(rankToIbis[i], bytesSentPerRank[i]);
		}
	}
	
	int getPort() {
		return serverSocket.socket().getLocalPort();
	}
	
	void close() {
		try {
			serverSocket.close();
		} catch (IOException e) {
			logger.error("Error while closing MPI profiling server socket", e);
		}
	}

	public void run() {
		while (serverSocket.isOpen()) {
			try {
				SocketChannel channel = serverSocket.accept();

				new MPIProfilingConnection(channel, this, rankToIbis.length);
			} catch (IOException e) {
				if (!serverSocket.isOpen()) {
					return;
				}
				logger.error("Error in accepting connection", e);
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e1) {
					return;
				}
			}
		}
	}
}
