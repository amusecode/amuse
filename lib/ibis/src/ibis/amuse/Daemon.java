package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCapabilities;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.PortType;
import ibis.ipl.RegistryEventHandler;

import java.io.IOException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Ibis Daemon. Starts workers for Amuse using software from the Ibis project.
 * 
 * @author Niels Drost
 * 
 */
public class Daemon implements RegistryEventHandler {

    private static class Shutdown extends Thread {
        private final Daemon daemon;

        Shutdown(Daemon daemon) {
            this.daemon = daemon;
        }

        public void run() {
            daemon.end();
        }
    }

    public static final int DEFAULT_PORT = 61575;

    public static PortType portType = new PortType(
            PortType.COMMUNICATION_RELIABLE, PortType.SERIALIZATION_OBJECT,
            PortType.RECEIVE_EXPLICIT, PortType.RECEIVE_TIMEOUT, PortType.CONNECTION_ONE_TO_ONE);

    public static IbisCapabilities ibisCapabilities = new IbisCapabilities(
            IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED,
            IbisCapabilities.TERMINATION);

    private static final Logger logger = LoggerFactory.getLogger(Daemon.class);

    private final ServerSocketChannel loopbackServer;

    private final Ibis ibis;

    private final Deployment deployment;

    public Daemon(int port) throws Exception {
        deployment = new Deployment();

        Properties properties = new Properties();
        properties.put("ibis.server.address", deployment.getServerAddress());
        properties.put("ibis.pool.name", "amuse");
        ibis = IbisFactory.createIbis(ibisCapabilities, properties, true, this,
                portType);

        loopbackServer = ServerSocketChannel.open();
        // bind to local host
        loopbackServer.socket().bind(
                new InetSocketAddress(InetAddress.getByName(null), port));

        ibis.registry().enableEvents();

        IbisIdentifier master = ibis.registry().elect("master");

        if (!ibis.identifier().equals(master)) {
            throw new IOException("did not win master election: another daemon"
                    + " must be present in this pool");
        }
    }

    public void run() {
        while (true) {
            SocketChannel socket = null;
            try {
                logger.info("Waiting for connection");
                socket = loopbackServer.accept();

                new LocalWorker(socket, ibis, deployment);
            } catch (Exception e) {
                if (socket != null) {
                    try {
                        socket.close();
                    } catch (Exception e2) {
                        // IGNORE
                    }
                }

                if (!loopbackServer.isOpen()) {
                    return;
                }
                logger.error("Error while handling new worker connection", e);

                // wait a bit before handling the next connection
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e1) {
                    // IGNORE
                }
            }
        }
    }

    public void end() {
        try {
            loopbackServer.close();
        } catch (IOException e) {
            logger.error("error on ending loopback server", e);
        }
        try {
            ibis.registry().terminate();

            ibis.end();
        } catch (IOException e) {
            logger.error("error on ending ibis", e);
        }
    }

    public static void main(String[] arguments) throws IOException {
        int port = DEFAULT_PORT;

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equals("-p") || arguments[i].equals("--port")) {
                i++;
                port = Integer.parseInt(arguments[i]);
            }
        }

        try {
            Daemon daemon = new Daemon(port);

            Runtime.getRuntime().addShutdownHook(new Shutdown(daemon));

            daemon.run();
        } catch (Exception e) {
            logger.error("Error while running daemon", e);
        }

    }

    @Override
    public void joined(IbisIdentifier joinedIbis) {
        logger.info("new ibis " + joinedIbis);
    }

    @Override
    public void left(IbisIdentifier leftIbis) {
        logger.info("ibis left " + leftIbis);
    }

    @Override
    public void died(IbisIdentifier corpse) {
        logger.info("ibis died " + corpse);
    }

    @Override
    public void gotSignal(String signal, IbisIdentifier source) {
        // IGNORE
    }

    @Override
    public void electionResult(String electionName, IbisIdentifier winner) {
        // IGNORE
    }

    @Override
    public void poolClosed() {
        // IGNORE
    }

    @Override
    public void poolTerminated(IbisIdentifier source) {
        logger.info("pool terminated!");
    }

}
