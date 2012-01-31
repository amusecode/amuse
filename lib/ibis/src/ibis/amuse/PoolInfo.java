package ibis.amuse;

import java.io.IOException;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCapabilities;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.IbisProperties;
import ibis.ipl.PortType;

/**
 * Keeps track of companion Ibisses of a parallel MPI code. Needed to initialize
 * MPI properly.
 * 
 * @author Niels Drost
 * 
 */
public class PoolInfo {

    private static class Shutdown extends Thread {
        private final Ibis ibis;

        Shutdown(Ibis ibis) {
            this.ibis = ibis;
        }

        public void run() {
            try {
                ibis.end();
            } catch (IOException e) {
                logger.error("Error ending PoolInfo");
            }
        }
    }

    private static final Logger logger = LoggerFactory.getLogger(PoolInfo.class);

    private final Ibis ibis;

    private final IbisIdentifier[] pool;
    private final String[] hostnames;
    private final int[] ports;
    private final int rank;
    private final int nrOfNodes;

    public static IbisCapabilities ibisCapabilities = new IbisCapabilities(IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED, IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS,
            IbisCapabilities.CLOSED_WORLD);

    PoolInfo(String poolID, int poolSize, int localMonitoringPort) throws Exception {
    	this.nrOfNodes = poolSize;
    	
        Properties properties = new Properties();
        properties.setProperty(IbisProperties.POOL_NAME, poolID);
        properties.setProperty(IbisProperties.POOL_SIZE, Integer.toString(poolSize));

        String host = InetAddress.getLocalHost().getHostAddress();

        ibis = IbisFactory.createIbis(ibisCapabilities, properties, true, null, null, host + "," + localMonitoringPort,
                new PortType[0]);

        // end ibis when shutdown
        Runtime.getRuntime().addShutdownHook(new Shutdown(ibis));

        logger.info("initializing poolinfo, waiting for others");

        ibis.registry().waitUntilPoolClosed();

        pool = ibis.registry().joinedIbises();

        hostnames = new String[pool.length];
        ports = new int[pool.length];

        int rank = 0;
        for (int i = 0; i < pool.length; i++) {
            String[] hostAndPort = pool[i].tagAsString().split(",");
            if (hostAndPort.length != 2) {
                throw new Exception("could not get host and port from string " + pool[i].tagAsString());
            }
            hostnames[i] = hostAndPort[0];
            ports[i] = Integer.parseInt(hostAndPort[1]);

            if (ibis.identifier().equals(pool[i])) {
                rank = i;
            }
        }
        this.rank = rank;

        logger.info("initialized poolinfo. Hosts: " + Arrays.toString(hostnames));
    }

    public String[] getHostnames() {
        return hostnames;
    }
    
    public int getNrOfNodes() {
    	return nrOfNodes;
    }

    public int[] getWorkerPortList(int nrOfWorkers) {
        int[] result = new int[nrOfWorkers];
        int next = 0;

        for (int i = 0; i < ports.length; i++) {
            // number of processes per node
            int workerCount = nrOfWorkers / ports.length;
            // nrOfWorkers not divideble by pool size. see if this is a
            // "remainder" node with an extra worker
            if (i < nrOfWorkers % ports.length) {
                workerCount++;
            }
            for (int j = 0; j < workerCount; j++) {
                result[next++] = ports[i];
            }
        }

        return result;
    }

    String[] getWorkerHostList(int nrOfWorkers) throws Exception {
        ArrayList<String> result = new ArrayList<String>();

        for (int i = 0; i < hostnames.length; i++) {
            // number of processes per node
            int workerCount = nrOfWorkers / hostnames.length;
            // nrOfWorkers not divideble by number of hosts. see if this is a
            // "remainder" node with an extra worker
            if (i < nrOfWorkers % hostnames.length) {
                workerCount++;
            }
            for (int j = 0; j < workerCount; j++) {
                result.add(hostnames[i]);
            }
        }
        if (result.size() != nrOfWorkers) {
            throw new Exception("error in creating host file. Size is " + result.size() + " but should be "
                    + nrOfWorkers);
        }

        return result.toArray(new String[0]);

    }

    IbisIdentifier[] getWorkerIbisList(int nrOfWorkers) throws Exception {
        ArrayList<IbisIdentifier> result = new ArrayList<IbisIdentifier>();

        for (int i = 0; i < pool.length; i++) {
            // number of processes per node
            int workerCount = nrOfWorkers / pool.length;
            // nrOfWorkers not divideble by number of hosts. see if this is a
            // "remainder" node with an extra worker
            if (i < nrOfWorkers % pool.length) {
                workerCount++;
            }
            for (int j = 0; j < workerCount; j++) {
                result.add(pool[i]);
            }
        }
        if (result.size() != nrOfWorkers) {
            throw new Exception("error in creating host file. Size is " + result.size() + " but should be "
                    + nrOfWorkers);
        }

        return result.toArray(new IbisIdentifier[0]);

    }

    public void terminate() {
        try {
            ibis.registry().terminate();
        } catch (IOException e) {
            logger.error("Could not terminate poolinfo pool", e);
        }
    }

    public void waitUntilTerminated() {
        ibis.registry().waitUntilTerminated();
    }

    public int getRank() {
        return rank;
    }

    public void end() throws IOException {
        ibis.end();
    }

}
