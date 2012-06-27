package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.util.RunProcess;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class representing a code that is used via a loopback socket interface.
 * 
 * @author Niels Drost
 * 
 */
public class SocketCodeInterface implements CodeInterface {

    private static final Logger logger = LoggerFactory.getLogger(SocketCodeInterface.class);

    private static final int ACCEPT_TRIES = 20;
    private static final int ACCEPT_TIMEOUT = 1000; // ms

    private static final String[] ENVIRONMENT_BLACKLIST = { "JOB_ID", "PE_", "PRUN_", "JOB_NAME", "JOB_SCRIPT", "OMPI_" };

    private File executable;

    // local socket communication stuff

    private ServerSocketChannel serverSocket;

    private SocketChannel socket;

    private Process process;
    private int result = 0;

    private AmuseMessage requestMessage;
    private AmuseMessage resultMessage;

    private OutputForwarder out;
    private OutputForwarder err;

    private final WorkerInfo info;
    private final AmuseConfigOptions amuseConfig;

    private final String[] hostnames;
    private final int[] ports;

    private final Ibis wideAreaIbis;

    SocketCodeInterface(WorkerInfo info, AmuseConfigOptions amuseConfig, IbisIdentifier[] localIbises, Ibis wideAreaIbis)
            throws CodeException {
        this.wideAreaIbis = wideAreaIbis;
        this.info = info;
        this.amuseConfig = amuseConfig;
        int next = 0;
        int nrOfProcesses = info.getNrOfProcesses();

        hostnames = new String[nrOfProcesses];
        ports = new int[info.getNrOfProcesses()];

        if (localIbises == null) {
            for (int i = 0; i < nrOfProcesses; i++) {
                hostnames[i] = "localhost";
                ports[i] = 0;
            }
        } else {
            for (int i = 0; i < localIbises.length; i++) {
                String[] hostAndPort = localIbises[i].tagAsString().split(",");
                if (hostAndPort.length != 2) {
                    throw new CodeException("could not get host and port from string " + localIbises[i].tagAsString());
                }
                String hostname = hostAndPort[0];
                int port = Integer.parseInt(hostAndPort[1]);

                // number of processes per node
                int workerCount = info.getNrOfProcesses() / localIbises.length;
                // nrOfWorkers not divideble by number of hosts. see if this is
                // a
                // "remainder" node with an extra worker
                if (i < info.getNrOfProcesses() % localIbises.length) {
                    workerCount++;
                }
                for (int j = 0; j < workerCount; j++) {
                    hostnames[next] = hostname;
                    ports[next] = port;

                    next++;
                }
            }
            if (next != info.getNrOfProcesses()) {
                logger.error("error in setting ibises. List is of length " + next + " but should be " + nrOfProcesses);
            }
        }
    }

    public void init() throws CodeException {
        if (info.copyWorkerCode()) {
            //executable in cwd
            executable = new File(info.getCodeName());
        } else {
            executable = new File(info.getAmuseHome() + File.separator + info.getCodeDir() + File.separator
                    + info.getCodeName());
        }
        if (!executable.isFile()) {
            throw new CodeException("Cannot find executable for code " + info.getCodeName() + ": " + executable);
        }

        if (!executable.canExecute() && ! info.copyWorkerCode()) {
            throw new CodeException(executable + " is not executable");
        }

        try {
            serverSocket = ServerSocketChannel.open();
            serverSocket.socket().bind(null);

            File hostFile = File.createTempFile("host", "file").getAbsoluteFile();

            FileWriter hostFileWriter = new FileWriter(hostFile);
            for (String hostname : hostnames) {
                hostFileWriter.write(hostname + "\n");
            }
            hostFileWriter.flush();
            hostFileWriter.close();

            logger.info("host file = " + hostFile);

            File portFile = File.createTempFile("port", "file").getAbsoluteFile();

            FileWriter portFileWriter = new FileWriter(portFile);
            for (int port : ports) {
                portFileWriter.write(port + "\n");
            }
            portFileWriter.flush();
            portFileWriter.close();

            logger.info("port file = " + portFile);

            ProcessBuilder builder = new ProcessBuilder();

            for (String key : builder.environment().keySet().toArray(new String[0])) {
                for (String blacklistedKey : ENVIRONMENT_BLACKLIST) {
                    if (key.startsWith(blacklistedKey)) {
                        builder.environment().remove(key);
                        logger.info("removed " + key + " from environment");
                    }
                }
            }

            builder.environment().put("OMPI_IBIS_PROFILING_PORT_FILE", portFile.getAbsolutePath());

            if (!amuseConfig.isMpiexecEnabled()) {
                logger.info("not using mpiexec (as it is disabled)");
                if (info.getNrOfProcesses() > 1) {
                    throw new CodeException(info.getNrOfProcesses() + "requested, but mpiexec disabled");
                }
            } else if (info.getNrOfNodes() == 1) {
                // no need for machine file, set number of processes.
                builder.command(amuseConfig.getMpiexec(), "-n", Integer.toString(info.getNrOfProcesses()));
            } else {
                // use machine file
                builder.command(amuseConfig.getMpiexec(), "-machinefile", hostFile.getAbsolutePath());
            }

            if (info.copyWorkerCode()) {
                // run executable via amuse.sh script to set python path and
                // interpreter
                builder.command().add(new File(amuseConfig.getAmuseHome(), "amuse.sh").getAbsolutePath());
            }

            // executable and port options
            builder.command().add(executable.toString());
            builder.command().add(Integer.toString(serverSocket.socket().getLocalPort()));

            // make sure there is an "output" directory for a code to put output
            // in
            new File("output").mkdir();

            logger.info("starting worker process, command = " + builder.command());

            synchronized (this) {
                process = builder.start();

                out = new OutputForwarder(process.getInputStream(), info.getStdoutFile(), wideAreaIbis);
                err = new OutputForwarder(process.getErrorStream(), info.getStderrFile(), wideAreaIbis);

            }

            logger.info("process started");

            socket = acceptConnection(serverSocket);

            logger.info("connection established");

        } catch (IOException e) {
            throw new CodeException("Cannot initialize socket code interface for " + info.getID() + ": "
                    + e.getMessage(), e);
        }
    }

    private static SocketChannel acceptConnection(ServerSocketChannel serverSocket) throws IOException {
        serverSocket.configureBlocking(false);
        for (int i = 0; i < ACCEPT_TRIES; i++) {
            SocketChannel result = serverSocket.accept();

            if (result != null) {
                return result;
            }
            try {
                Thread.sleep(ACCEPT_TIMEOUT);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }
        throw new IOException("worker not started, socket connection failed to initialize");
    }

    @Override
    public void call() throws CodeException {
        try {
            logger.debug("performing call with function ID " + requestMessage.getFunctionID());
            requestMessage.writeTo(socket);
            resultMessage.readFrom(socket);
            logger.debug("done performing call with function ID " + requestMessage.getFunctionID() + " error = "
                    + resultMessage.getError());
        } catch (IOException e) {
            throw new CodeException("communication with code failed", e);
        }
    }

    @Override
    public synchronized void end() {
        if (process != null) {
            process.destroy();

            try {
                result = process.exitValue();
                logger.info("Process ended with result " + result);
            } catch (IllegalThreadStateException e) {
                logger.error("Process not ended after process.destroy()!");
            }
        }

        if (out != null) {
            // wait for out and err a bit
            try {
                out.join(1000);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }

        if (err != null) {
            try {
                err.join(1000);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }
    }

    @Override
    public void setRequestMessage(AmuseMessage message) {
        this.requestMessage = message;
    }

    @Override
    public void setResultMessage(AmuseMessage message) {
        this.resultMessage = message;
    }

    @Override
    public synchronized int getResult() {
        return result;
    }

}
