package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCapabilities;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.IbisProperties;
import ibis.ipl.PortType;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.net.InetAddress;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Proxy to a code running on some resource. Handles all communication with the
 * Worker running on the machine running amuse. Talks to the actual code using a
 * CodeInterface
 * 
 * @author Niels Drost
 * 
 */
public class CodeProxy {

    private static class Shutdown extends Thread {
        private final CodeProxy codeProxy;

        Shutdown(CodeProxy codeProxy) {
            this.codeProxy = codeProxy;
        }

        public void run() {
            codeProxy.end();
        }
    }

    private static Logger logger = LoggerFactory.getLogger(CodeProxy.class);

    public static IbisCapabilities localIbisCapabilities = new IbisCapabilities(IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED, IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS,
            IbisCapabilities.CLOSED_WORLD);

    private final AmuseConfigOptions amuseConfig;

    // collects mpi profiling data. only used if multiple nodes are used.
    private final MPIProfilingCollector mpiProfilingCollector;

    // ibis used to communicate with other nodes running this code.
    private final Ibis localIbis;

    // list of all ibises running this code.
    private final IbisIdentifier[] localIbises;

    private final boolean master;

    private final SendPort sendPort;

    private final ReceivePort receivePort;

    private final Ibis wideAreaIbis;

    private final CodeInterface codeInterface;

    private boolean ended = false;

    public CodeProxy(WorkerInfo info) throws Exception {
        this.amuseConfig = new AmuseConfigOptions(info.getAmuseHome());

        // hostname used for MPI
        String hostname = InetAddress.getLocalHost().getHostAddress();

        if (info.getNrOfNodes() == 1) {
            // no need to measure mpi traffic
            mpiProfilingCollector = null;

            // No need for a local Ibis, there is only one node (us)
            localIbis = null;
            localIbises = null;
            master = true;
        } else {
            // start monitoring if mpi profiling data
            mpiProfilingCollector = new MPIProfilingCollector();

            // create Ibis, wait for everyone to join, create list of
            // identifiers
            Properties properties = new Properties();
            properties.put(IbisProperties.POOL_NAME, info.getID());
            properties.put(IbisProperties.POOL_SIZE, Integer.toString(info.getNrOfNodes()));

            localIbis = IbisFactory.createIbis(localIbisCapabilities, properties, true, null, null, hostname + ","
                    + mpiProfilingCollector.getPort(), new PortType[0]);
            localIbis.registry().waitUntilPoolClosed();
            localIbises = localIbis.registry().joinedIbises();

            // we are master if we are first
            master = (localIbis.identifier().equals(localIbises[0]));
            mpiProfilingCollector.setIbises(localIbises, info.getNrOfProcesses());
        }

        if (master) {
            logger.info("Master = " + hostname);
            wideAreaIbis = IbisFactory.createIbis(Daemon.ibisCapabilities, new Properties(), true, null,
                    Daemon.portType);
            sendPort = wideAreaIbis.createSendPort(Daemon.portType);
            receivePort = wideAreaIbis.createReceivePort(Daemon.portType, info.getID());

            receivePort.enableConnections();
            //wideAreaIbis.registry().enableEvents();

            // get address of amuse node
            logger.debug("waiting for result of amuse election");
            IbisIdentifier amuse = wideAreaIbis.registry().getElectionResult("amuse");

            logger.debug("connecting to receive port of worker at amuse node");
            sendPort.connect(amuse, info.getID());
            logger.debug("connected, saying hello");
            WriteMessage message = sendPort.newMessage();
            message.writeObject(receivePort.identifier());
            message.finish();

            if (info.getInterfaceType() == null || info.getInterfaceType().equals("sockets")) {
                codeInterface = new SocketCodeInterface(info, amuseConfig, localIbises, wideAreaIbis);
            } else if (info.getInterfaceType().equals("jni")) {
                codeInterface = new JNICodeInterface(info);
            } else {
                throw new CodeException("unknown worker interface type: " + info.getInterfaceType());
            }
        } else {
            wideAreaIbis = null;
            sendPort = null;
            receivePort = null;
            codeInterface = null;
        }
    }

    private boolean isMaster() {
        return master;
    }
    

    private int getResult() {
        return codeInterface.getResult();
    }

    protected AmuseMessage receiveInitRequest() throws IOException {
        AmuseMessage result = new AmuseMessage();

        ReadMessage message = receivePort.receive();
        result.readFrom(message);
        message.finish();

        if (result.getFunctionID() != AmuseMessage.FUNCTION_ID_INIT) {
            throw new IOException("first call must be init function");
        }

        return result;
    }

    protected void sendInitReply(int callID) throws IOException {
        AmuseMessage result = new AmuseMessage(callID, AmuseMessage.FUNCTION_ID_INIT, 1);

        WriteMessage message = sendPort.newMessage();
        result.writeTo(message);
        message.finish();
    }

    protected void sendInitReply(int callID, Throwable error) throws IOException {
        AmuseMessage result = new AmuseMessage(callID, AmuseMessage.FUNCTION_ID_INIT, 1, error.getMessage());

        WriteMessage message = sendPort.newMessage();
        result.writeTo(message);
        message.finish();
    }

    synchronized void end() {
        if (ended) {
            return;
        }
        ended = true;

        if (codeInterface != null) {
            codeInterface.end();
        }

        if (sendPort != null) {
            try {
                sendPort.close();
            } catch (IOException e) {
                logger.error("Error on ending code interface", e);
            }
        }

        if (receivePort != null) {
            try {
                receivePort.close(1000);
            } catch (IOException e) {
                logger.error("Error on ending code interface", e);
            }
        }

        if (wideAreaIbis != null) {
            try {
                wideAreaIbis.end();
            } catch (IOException e) {
                logger.error("Error on ending code interface", e);
            }
        }

        if (localIbis != null) {
            if (master) {
                try {
                    localIbis.registry().terminate();
                } catch (IOException e) {
                    logger.error("Error on ending code interface", e);
                }
            }
            try {
                localIbis.end();
            } catch (IOException e) {
                logger.error("Error on ending code interface", e);
            }
        }
    }

    /**
     * Continuously receives a message, performs a call, sends a reply.
     */
    public void runMaster() {
        boolean running = true;
        long start, finish;
        AmuseMessage requestMessage = new AmuseMessage();
        AmuseMessage resultMessage = new AmuseMessage();

        while (running) {
            start = System.currentTimeMillis();
            try {
                logger.debug("Receiving call message");
                ReadMessage readMessage = receivePort.receive();

                if (readMessage == null) {
                    throw new IOException("cannot get request from worker");
                }

                logger.debug("Reading call request from IPL message");

                boolean changed = requestMessage.readFrom(readMessage);

                readMessage.finish();

                if (changed) {
                    // read method indicates one or more buffers has changed in
                    // the message
                    logger.debug("(re)setting request message");
                    codeInterface.setRequestMessage(requestMessage);
                }

                int functionID = requestMessage.getFunctionID();

                if (functionID == AmuseMessage.FUNCTION_ID_STOP) {
                    // final request handled
                    running = false;
                }

                // clean message
                resultMessage.clear();

                logger.debug("Performing call for function " + functionID);
                try {
                    if (functionID == AmuseMessage.FUNCTION_ID_INIT) {
                        codeInterface.setRequestMessage(requestMessage);
                        codeInterface.setResultMessage(resultMessage);
                        codeInterface.init();
                    } else {
                        // perform call. Will put result in result message
                        codeInterface.call();
                    }
                } catch (CodeException e) {
                    logger.error("exception on performing call", e);
                    // put an exception in the result message
                    resultMessage.clear();
                    resultMessage.setCallID(requestMessage.getCallID());
                    resultMessage.setFunctionID(requestMessage.getFunctionID());
                    resultMessage.setCallCount(requestMessage.getCount());
                    resultMessage.setError(e.getMessage());
                }

                logger.debug("result: " + resultMessage);

                WriteMessage writeMessage = sendPort.newMessage();

                resultMessage.writeTo(writeMessage);

                writeMessage.finish();

                logger.debug("Done performing call for function " + functionID);
                finish = System.currentTimeMillis();

                if (logger.isDebugEnabled()) {
                    logger.debug("Call took " + (finish - start) + " ms");
                }
            } catch (Exception e) {
                logger.error("Error while handling request", e);
                running = false;
            }
        }
    }

    public void runSlave() {
        // Wait until the local Ibis is terminated, then exit
        localIbis.registry().waitUntilTerminated();
    }

    public static void main(String[] arguments) {
        WorkerInfo info = new WorkerInfo();

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equals("--worker-id")) {
                i++;
                info.setWorkerID(arguments[i]);
            } else if (arguments[i].equals("--code-name")) {
                i++;
                info.setCodeName(arguments[i]);
            } else if (arguments[i].equals("--code-dir")) {
                i++;
                info.setCodeDir(arguments[i]);
            } else if (arguments[i].equals("--amuse-home")) {
                i++;
                info.setAmuseHome(arguments[i]);
            } else if (arguments[i].equals("--interface-type")) {
                i++;
                info.setInterfaceType(arguments[i]);
            } else if (arguments[i].equals("--number-of-processes")) {
                i++;
                info.setNrOfProcesses(Integer.parseInt(arguments[i]));
            } else if (arguments[i].equals("--number-of-nodes")) {
                i++;
                info.setNrOfNodes(Integer.parseInt(arguments[i]));
            } else if (arguments[i].equals("--number-of-threads")) {
                i++;
                info.setNrOfThreads(Integer.parseInt(arguments[i]));
            } else if (arguments[i].equals("--stdout")) {
                i++;
                info.setStdoutFile(arguments[i]);
            } else if (arguments[i].equals("--stderr")) {
                i++;
                info.setStderrFile(arguments[i]);
            } else if (arguments[i].equals("--python-code")) {
                i++;
                info.setIsPythonCode(Boolean.parseBoolean(arguments[i]));
            } else if (arguments[i].equals("--copy-worker-code")) {
                i++;
                info.setCopyWorkerCode(Boolean.parseBoolean(arguments[i]));
            } else {
                System.err.println("Unknown option: " + arguments[i]);
                System.exit(1);
            }
        }

        try {
            CodeProxy codeProxy = new CodeProxy(info);

            Runtime.getRuntime().addShutdownHook(new Shutdown(codeProxy));
            
            if (codeProxy.isMaster()) {
                codeProxy.runMaster();
            } else {
                codeProxy.runSlave();
            }
            codeProxy.end();
            System.exit(codeProxy.getResult());
        } catch (Exception e) {
            logger.error("Error on running code proxy", e);
            System.exit(1);
        }
    }
}
