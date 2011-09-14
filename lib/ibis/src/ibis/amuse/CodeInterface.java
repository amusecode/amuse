package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.RegistryEventHandler;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class CodeInterface implements RegistryEventHandler, Runnable {

    private static class Shutdown extends Thread {
        private final CodeInterface codeInterface;

        Shutdown(CodeInterface worker) {
            this.codeInterface = worker;
        }

        public void run() {
            codeInterface.end();
        }
    }

    private static Logger logger = LoggerFactory.getLogger(CodeInterface.class);

    private final PoolInfo poolInfo;

    private final SendPort sendPort;

    private final ReceivePort receivePort;

    private final Ibis ibis;

    private volatile boolean ended = false;

    private final AmuseMessage requestMessage;
    private final AmuseMessage resultMessage;

    public CodeInterface(String id, PoolInfo poolInfo) throws Exception {

        this.poolInfo = poolInfo;

        requestMessage = new AmuseMessage();
        resultMessage = new AmuseMessage();

        ibis = IbisFactory.createIbis(Daemon.ibisCapabilities, this,
                Daemon.portType);

        sendPort = ibis.createSendPort(Daemon.portType);
        receivePort = ibis.createReceivePort(Daemon.portType, id);

        try {
            receivePort.enableConnections();
            ibis.registry().enableEvents();

            // get address of amuse node
            logger.debug("waiting for result of amuse election");
            IbisIdentifier amuse = ibis.registry().getElectionResult("amuse");

            logger.debug("connecting to receive port of worker at amuse node");
            sendPort.connect(amuse, id);
            logger.debug("connected, saying hello");
            WriteMessage message = sendPort.newMessage();
            message.writeObject(receivePort.identifier());
            message.finish();

        } catch (Exception e) {
            end();
            throw e;
        }
        logger.info("Generic Code interface initialized");
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

    void end() {
        if (ended) {
            return;
        }
        ended = true;

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

        if (ibis != null) {
            try {
                ibis.end();
            } catch (IOException e) {
                logger.error("Error on ending code interface", e);
            }
        }

        if (poolInfo != null) {
            poolInfo.terminate();
        }
    }

    @Override
    public void joined(IbisIdentifier joinedIbis) {
        // IGNORE
    }

    @Override
    public void left(IbisIdentifier leftIbis) {
        // IGNORE
    }

    @Override
    public void died(IbisIdentifier corpse) {
        // IGNORE
    }

    @Override
    // signal send by amuse/ibis daemon to THIS worker
    public void gotSignal(String signal, IbisIdentifier source) {
        if (signal.equals("end")) {
            logger.info("got end signal");
            end();
        }
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
    // signal send by amuse/ibis daemon to ALL workers
    public void poolTerminated(IbisIdentifier source) {
        end();
    }

    /**
     * Continuously receives a message, performs a call, sends a reply.
     */
    public void run() {
        boolean running = true;
        long start, finish;
        
        setRequestMessage(requestMessage);
        setResultMessage(resultMessage);

        while (running) {
            start = System.currentTimeMillis();
            try {
                logger.debug("Receiving call message");
                ReadMessage readMessage = receivePort.receive();

                logger.debug("Reading call request from IPL message");

                boolean changed = requestMessage.readFrom(readMessage);

                readMessage.finish();

                if (changed) {
                    // read method indicates one or more buffers has changed in
                    // the message
                    logger.debug("(re)setting request message");
                    setRequestMessage(requestMessage);
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
                    // perform call. Will put result in result message
                    call();
                } catch (Exception e) {
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
            } catch (Throwable e) {
                logger.error("Error while handling request", e);
                try {
                    Thread.sleep(1000);
                } catch (Exception e2) {
                    // IGNORE
                }
            }
        }

    }

    /**
     * Sets the message used as input in a call. Also used when the buffers
     * inside the message change
     * 
     * @param message
     *            the message
     */
    abstract void setRequestMessage(AmuseMessage message);

    /**
     * Sets the message used as output in a call
     * 
     * @param message
     *            the message
     */
    abstract void setResultMessage(AmuseMessage message);

    /**
     * Performs a call. Uses data from the request and result message
     * 
     * @throws Exception
     *             if the function does not exist, or any other error occurs.
     */
    abstract void call() throws Exception;
    
    public static void main(String[] arguments) {
        String id = null;
        String codeName = null;
        boolean jni = false;
        String amuseHome = null;
        String codeDir = null;
        int nrOfWorkers = 1;
        int nrOfNodes = 1;
        String mpirun = null;

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equals("-i") || arguments[i].equals("--worker-id")) {
                i++;
                id = arguments[i];
            } else if (arguments[i].equals("-c")
                    || arguments[i].equals("--code-name")) {
                i++;
                codeName = arguments[i];
            } else if (arguments[i].equals("-d")
                    || arguments[i].equals("--code-dir")) {
                i++;
                codeDir = arguments[i];
            } else if (arguments[i].equals("-a")
                    || arguments[i].equals("--amuse-home")) {
                i++;
                amuseHome = arguments[i];
            } else if (arguments[i].equals("-j")
                    || arguments[i].equals("--jni")) {
                jni = true;
            } else if (arguments[i].equals("-n")
                    || arguments[i].equals("--number-of-workers")) {
                i++;
                nrOfWorkers = Integer.parseInt(arguments[i]);
            } else if (arguments[i].equals("--number-of-nodes")) {
                i++;
                nrOfNodes = Integer.parseInt(arguments[i]);
            } else if (arguments[i].equals("-m")
                    || arguments[i].equals("--mpirun")) {
                i++;
                mpirun = arguments[i];
            } else {
                System.err.println("Unknown option: " + arguments[i]);
                System.exit(1);
            }
        }

        CodeInterface codeInterface = null;
        PoolInfo poolInfo = null;

        try {
            logger.info("Starting worker " + id + " running " + codeName);

            poolInfo = new PoolInfo(id, nrOfNodes);

            if (poolInfo.getRank() == 0) {
                if (jni) {
                    codeInterface = new JNICodeInterface(id, poolInfo, codeName);
                } else {
                    codeInterface = new SocketCodeInterface(id, poolInfo, codeName,
                            codeDir, amuseHome, mpirun, nrOfWorkers);
                }

                Runtime.getRuntime().addShutdownHook(
                        new Shutdown(codeInterface));

                codeInterface.run();

                logger.info("Worker " + id + " running " + codeName + " DONE!");

                poolInfo.terminate();
            } else {
                poolInfo.waitUntilTerminated();
            }
            poolInfo.end();
        } catch (Exception e) {
            logger.error("Error on running code server", e);
        } finally {
            if (codeInterface != null) {
                codeInterface.end();
            }
            if (poolInfo != null) {
                poolInfo.terminate();
            }
        }
    }

}
