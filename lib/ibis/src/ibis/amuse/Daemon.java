package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCapabilities;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.PortType;
import ibis.ipl.RegistryEventHandler;

import java.io.File;
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

    public static PortType portType = new PortType(PortType.COMMUNICATION_RELIABLE, PortType.SERIALIZATION_OBJECT,
            PortType.RECEIVE_EXPLICIT, PortType.RECEIVE_TIMEOUT, PortType.CONNECTION_ONE_TO_ONE);

    public static IbisCapabilities ibisCapabilities = new IbisCapabilities(IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED, IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS);

    private static final Logger logger = LoggerFactory.getLogger(Daemon.class);

    private final ServerSocketChannel loopbackServer;

    private final Ibis ibis;

    private final Deployment deployment;

    public Daemon(int port, boolean verbose, boolean keepSandboxes, boolean gui, boolean hubs, File jungle)
            throws Exception {
        deployment = new Deployment(verbose, keepSandboxes, gui, hubs, jungle);

        Properties properties = new Properties();
        properties.put("ibis.server.address", deployment.getServerAddress());
        properties.put("ibis.pool.name", "amuse");
        properties.put("ibis.location", "daemon@local");
        properties.put("ibis.managementclient", "true");
        properties.put("ibis.bytescount", "true");

        ibis = IbisFactory.createIbis(ibisCapabilities, properties, true, this, portType);

        loopbackServer = ServerSocketChannel.open();
        // bind to local host
        loopbackServer.socket().bind(new InetSocketAddress(InetAddress.getByName(null), port));

        ibis.registry().enableEvents();

        IbisIdentifier amuse = ibis.registry().elect("amuse");

        if (!ibis.identifier().equals(amuse)) {
            throw new IOException("did not win amuse election: another daemon" + " must be present in this pool");
        }
        logger.info("Daemon running on port " + port);
    }

    public void run() {
        while (true) {
            SocketChannel socket = null;
            try {
                logger.debug("Waiting for connection");
                socket = loopbackServer.accept();

                new RemoteCodeInterface(socket, ibis, deployment);
            } catch (Exception e) {
                // report error to amuse
                AmuseMessage errormessage = new AmuseMessage(0, AmuseMessage.FUNCTION_ID_INIT, 1, e.getMessage());
                try {
                    errormessage.writeTo(socket);
                } catch (Exception e2) {
                    // IGNORE
                }
                logger.error("error on starting remote code", e);
                
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

    public static void printUsage() {
        System.err.println("Usage: ibis-deploy.sh [OPTIONS]\n" + "Options:\n"
                + "-p | --port [PORT]\t\tPort to listen on (default: " + DEFAULT_PORT + ")\n"
                + "-v | --verbose\t\t\tBe more verbose\n" + "-g | --gui\t\t\tStart a monitoring gui as well\n"
                + "-k | --keep-sandboxes\t\tKeep JavaGAT sandboxes (mostly for debugging)\n"
                + "-j | --jungle-file [FILE]\tName of jungle configuration file\n" + "-h | --help\t\t\tThis message");
    }

    public static void main(String[] arguments) throws IOException {
        int port = DEFAULT_PORT;
        boolean verbose = false;
        boolean gui = false;
        boolean hubs = true;
        boolean keepSandboxes = false;
        File jungle = null;

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equals("-p") || arguments[i].equals("--port")) {
                i++;
                port = Integer.parseInt(arguments[i]);
            } else if (arguments[i].equals("-v") || arguments[i].equals("--verbose")) {
                verbose = true;
            } else if (arguments[i].equals("-g") || arguments[i].equals("--gui")) {
                gui = true;
            } else if (arguments[i].equals("-k") || arguments[i].equals("--keep-sandboxes")) {
                keepSandboxes = true;
            } else if (arguments[i].equals("--no-hubs")) {
                hubs = false;
            } else if (arguments[i].equals("-j") || arguments[i].equals("--jungle-file")) {
                i++;
                jungle = new File(arguments[i]);
            } else if (arguments[i].equals("-h") || arguments[i].equals("--help")) {
                printUsage();
                return;
            } else {
                System.err.println("Unknown option: " + arguments[i]);
                printUsage();
                System.exit(1);
            }
        }

        try {
            Daemon daemon = new Daemon(port, verbose, keepSandboxes, gui, hubs, jungle);

            Runtime.getRuntime().addShutdownHook(new Shutdown(daemon));

            daemon.run();
        } catch (Exception e) {
            logger.error("Error while running daemon", e);
        }

    }

    @Override
    public void joined(IbisIdentifier joinedIbis) {
        logger.debug("new ibis " + joinedIbis);
    }

    @Override
    public void left(IbisIdentifier leftIbis) {
        logger.debug("ibis left " + leftIbis);
    }

    @Override
    public void died(IbisIdentifier corpse) {
        logger.error("ibis died " + corpse);
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
