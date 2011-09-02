package ibis.amuse;

import java.io.IOException;
import java.net.InetAddress;
import java.util.Arrays;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCapabilities;
import ibis.ipl.IbisCreationFailedException;
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
    
    private static final Logger logger = LoggerFactory.getLogger(PoolInfo.class);

    private final Ibis ibis;

    private final IbisIdentifier[] pool;
    private final String[] hostnames;
    private final int rank;

    public static IbisCapabilities ibisCapabilities = new IbisCapabilities(
            IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED,
            IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS,
            IbisCapabilities.CLOSED_WORLD);

    PoolInfo(String poolID, int poolSize) throws Exception {
        Properties properties = new Properties();

        properties.setProperty(IbisProperties.POOL_NAME, poolID);
        properties.setProperty(IbisProperties.POOL_SIZE, Integer.toString(poolSize));
        
        String hostname = InetAddress.getLocalHost().getHostName();
        
        ibis = IbisFactory.createIbis(ibisCapabilities, properties, true, null, null, hostname, new PortType[0]);
        
        logger.info("initialized poolinfo, waiting for others");

        ibis.registry().waitUntilPoolClosed();

        pool = ibis.registry().joinedIbises();

        hostnames = new String[pool.length];

        int rank = 0;
        for (int i = 0; i < hostnames.length; i++) {
            hostnames[i] = pool[i].tagAsString();

            if (ibis.identifier().equals(pool[i])) {
                rank = i;
            }
        }
        this.rank = rank;
        
        System.err.println("Hostnames: " + Arrays.toString(hostnames));
    }

    public String[] getHostnames() {
        return hostnames;
    }

    public void terminate() throws IOException {
        ibis.registry().terminate();
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
