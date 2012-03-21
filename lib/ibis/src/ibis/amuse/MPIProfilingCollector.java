package ibis.amuse;

import ibis.ipl.IbisIdentifier;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;
import java.util.HashMap;
import java.util.Map;

import javax.management.MBeanServer;
import javax.management.ObjectName;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MPIProfilingCollector extends Thread implements MPIProfilingCollectorMBean {

    private static final Logger logger = LoggerFactory.getLogger(MPIProfilingCollector.class);

    private final Map<IbisIdentifier, Long> sent;

    private IbisIdentifier[] rankToIbis = null;

    // connection with codes
    private final ServerSocketChannel serverSocket;

    public MPIProfilingCollector() throws Exception {
        sent = new HashMap<IbisIdentifier, Long>();

        MBeanServer mbs = ManagementFactory.getPlatformMBeanServer();
        ObjectName name = new ObjectName("ibis.amuse:type=MPIProfilingCollector");
        mbs.registerMBean(this, name);

        serverSocket = ServerSocketChannel.open();
        serverSocket.socket().bind(null);

        this.setDaemon(true);
        this.setName(this.getClass().getName());
        this.start();
    }

    synchronized void setIbises(IbisIdentifier[] ibises, int nrOfProcesses) {
        rankToIbis = new IbisIdentifier[nrOfProcesses];
        int next = 0;

        for (int i = 0; i < ibises.length; i++) {
            // number of processes per node
            int workerCount = nrOfProcesses / ibises.length;
            // nrOfWorkers not divideble by number of hosts. see if this is a
            // "remainder" node with an extra worker
            if (i < nrOfProcesses % ibises.length) {
                workerCount++;
            }
            for (int j = 0; j < workerCount; j++) {
                rankToIbis[next] = ibises[i];
                next++;
            }
        }
        if (next != rankToIbis.length) {
            logger.error("error in setting ibises. List is of length " + next + " but should be "
                    + nrOfProcesses);
        }
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
        if (rankToIbis == null) {
            return;
        }
        for (int i = 0; i < bytesSentPerRank.length; i++) {
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

                if (rankToIbis == null) {
                    throw new IOException("unknown pool size");
                }

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
