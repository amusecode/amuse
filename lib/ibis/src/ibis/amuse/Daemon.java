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
import java.nio.ByteBuffer;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;
import java.util.ArrayList;
import java.util.Properties;

import org.apache.log4j.FileAppender;
import org.apache.log4j.HTMLLayout;
import org.apache.log4j.PatternLayout;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Ibis Daemon. Starts workers for Amuse using software from the Ibis project.
 * 
 * @author Niels Drost
 * 
 */
public class Daemon implements RegistryEventHandler {

    // magic strings to denote the type of connection (and to make sure we are
    // talking to amuse)

    public static final String WORKER_TYPE_STRING = "TYPE_WORKER";

    public static final String OUTPUT_TYPE_STRING = "TYPE_OUTPUT";

    public static final int TYPE_STRING_LENGTH = WORKER_TYPE_STRING.length();

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

    public static final File DEFAULT_LOG_DIR = new File(System.getProperty("java.io.tmpdir"), "ibis-amuse-logs");

    public static PortType portType = new PortType(PortType.COMMUNICATION_RELIABLE, PortType.SERIALIZATION_OBJECT,
            PortType.RECEIVE_EXPLICIT, PortType.RECEIVE_TIMEOUT, PortType.CONNECTION_MANY_TO_ONE);

    public static IbisCapabilities ibisCapabilities = new IbisCapabilities(IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED, IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS);

    private static final Logger logger = LoggerFactory.getLogger(Daemon.class);

    private final ServerSocketChannel loopbackServer;

    private final Ibis ibis;

    private final Deployment deployment;
    
    private final OutputManager outputManager;

    public Daemon(int port, boolean verbose, boolean keepSandboxes, boolean gui, boolean useHubs, File[] jungleFiles,
            String[] hubs, File logDir) throws Exception {
        logDir.mkdirs();
        if (!logDir.isDirectory()) {
            throw new Exception("log dir (" + logDir + ") is not a directory, and cannot be created");
        }

        // use log4j for logger configuration: write logs to html file
        org.apache.log4j.Logger.getRootLogger().addAppender(
                new FileAppender(new HTMLLayout(), new File(logDir, "deploy.log.html").getPath(), false));

        // use log4j for logger configuration: write logs to txt file
        org.apache.log4j.Logger.getRootLogger().addAppender(
                new FileAppender(new PatternLayout("%d{HH:mm:ss} %-5p [%t] %c - %m%n"), new File(logDir,
                        "deploy.log.txt").getPath(), false));

        deployment = new Deployment(verbose, keepSandboxes, gui, useHubs, jungleFiles, hubs, logDir);

        Properties properties = new Properties();
        properties.put("ibis.server.address", deployment.getServerAddress());
        properties.put("ibis.pool.name", "amuse");
        properties.put("ibis.location", "daemon@local");
        properties.put("ibis.managementclient", "true");
        properties.put("ibis.bytescount", "true");

        ibis = IbisFactory.createIbis(ibisCapabilities, properties, true, this, portType);
        
        outputManager = new OutputManager(ibis);

        loopbackServer = ServerSocketChannel.open();
        // bind to local host
        loopbackServer.socket().bind(new InetSocketAddress(InetAddress.getByName(null), port));

        ibis.registry().enableEvents();

        IbisIdentifier amuse = ibis.registry().elect("amuse");

        if (!ibis.identifier().equals(amuse)) {
            throw new IOException("did not win amuse election: another daemon" + " must be present in this pool");
        }
        logger.info("Daemon running on port " + port);
        logger.info("Logging worker output to " + logDir);
    }

    public void run() {
        while (true) {
            SocketChannel socket = null;
            try {
                logger.debug("Waiting for connection");
                socket = loopbackServer.accept();
                logger.debug("New connection accepted");

                // read string, to make sure we are talking to amuse, and to get
                // the type of connection

                ByteBuffer magic = ByteBuffer.allocate(TYPE_STRING_LENGTH);

                while (magic.hasRemaining()) {
                    long read = socket.read(magic);

                    if (read == -1) {
                        throw new IOException("Connection closed on reading magic string");
                    }
                }

                String receivedString = new String(magic.array(), "UTF-8");
                if (receivedString.equalsIgnoreCase(WORKER_TYPE_STRING)) {
                    new Worker(socket, ibis, deployment);
                } else if (receivedString.equalsIgnoreCase(OUTPUT_TYPE_STRING)) {                
                    outputManager.newOutputConnection(socket);
                } else {
                    throw new IOException("magic string (" + WORKER_TYPE_STRING + " or " + OUTPUT_TYPE_STRING + ") not received. Instead got: "
                            + receivedString);
                }

                logger.debug("New connection handled");
            } catch (Exception e) {
                if (socket != null) {
                    // report error to amuse
                    AmuseMessage errormessage = new AmuseMessage(0, AmuseMessage.FUNCTION_ID_INIT, 1, e.getMessage());
                    try {
                        errormessage.writeTo(socket);
                    } catch (Exception e2) {
                        // IGNORE
                    }
                    logger.error("error on starting remote code", e);

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
                + "-j | --jungle-file [FILE]\tName of jungle configuration file\n"
                + "-h | --hub [RESOURCE_NAME]\tStart a hub on the given resource\n"
                + "-l | --log-dir [DIR]\t\tSet location of logging files\n" + "\t\t\t\t(defaults to \""
                + DEFAULT_LOG_DIR + "\")\n" + "-n | --no-hubs\t\t\tDo not start any hubs automatically\n"
                + "\t\t\t\t(except when specified via --hub option)\n" + "-h | --help\t\t\tThis message");
    }

    public static void main(String[] arguments) throws IOException {
        int port = DEFAULT_PORT;
        boolean verbose = false;
        boolean gui = false;
        boolean startHubs = true;
        boolean keepSandboxes = false;
        ArrayList<File> jungleFiles = new ArrayList<File>();
        ArrayList<String> hubs = new ArrayList<String>();
        File logDir = DEFAULT_LOG_DIR;

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
            } else if (arguments[i].equals("-n") || arguments[i].equals("--no-hubs")) {
                startHubs = false;
            } else if (arguments[i].equals("-j") || arguments[i].equals("--jungle-file")) {
                i++;
                jungleFiles.add(new File(arguments[i]));
            } else if (arguments[i].equals("-h") || arguments[i].equals("--hub")) {
                i++;
                hubs.add(arguments[i]);
            } else if (arguments[i].equals("-l") || arguments[i].equals("--log-dir")) {
                i++;
                logDir = new File(arguments[i]);
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
            Daemon daemon = new Daemon(port, verbose, keepSandboxes, gui, startHubs, jungleFiles.toArray(new File[0]),
                    hubs.toArray(new String[0]), logDir);

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
