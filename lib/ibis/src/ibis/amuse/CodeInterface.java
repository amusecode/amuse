package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.IbisProperties;
import ibis.ipl.ReceivePort;
import ibis.ipl.RegistryEventHandler;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;
import ibis.util.TypedProperties;

import java.io.IOException;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CodeInterface implements RegistryEventHandler {

    private static class Shutdown extends Thread {
        private final CodeInterface worker;

        Shutdown(CodeInterface worker) {
            this.worker = worker;
        }

        public void run() {
            worker.end();
        }
    }

    private static Logger logger = LoggerFactory.getLogger(CodeInterface.class);

    private final Runnable communityCode;

    private final SendPort sendPort;

    private final ReceivePort receivePort;

    private final Ibis ibis;

    private volatile boolean ended = false;

    public CodeInterface(String id, String codeName, String codeDir,
            String amuseHome, boolean jni, String[] hostnames) throws Exception {

        ibis = IbisFactory.createIbis(Daemon.ibisCapabilities, this,
                Daemon.portType);

        sendPort = ibis.createSendPort(Daemon.portType);
        receivePort = ibis.createReceivePort(Daemon.portType, id);

        try {

            if (jni) {
                communityCode = new JNICodeInterface(codeName, receivePort,
                        sendPort);
            } else {
                communityCode = new SocketCodeInterface(codeName, codeDir,
                        amuseHome, receivePort, sendPort);
            }

            receivePort.enableConnections();
            ibis.registry().enableEvents();

            // get address of amuse node
            logger.debug("waiting for result of amuse election");
            IbisIdentifier amuse = ibis.registry().getElectionResult("amuse");

            // connect to receive port of worker at amuse node, send hello
            sendPort.connect(amuse, id);
            WriteMessage message = sendPort.newMessage();
            message.writeObject(receivePort.identifier());
            message.finish();

        } catch (Exception e) {
            end();
            throw e;
        }
    }

    public void end() {
        if (ended) {
            return;
        }
        ended = true;

        try {
            sendPort.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        try {
            receivePort.close(1000);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        try {
            ibis.end();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
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
    public void poolTerminated(IbisIdentifier source) {
        end();
    }

    public void run() {
        communityCode.run();
        end();
    }

    public static void main(String[] arguments) {
        String id = null;
        String codeName = null;
        boolean jni = false;
        String amuseHome = null;
        String codeDir = null;
        int nrOfWorkers = 1;

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
            } else {
                System.err.println("Unknown option: " + arguments[i]);
                System.exit(1);
            }
        }

        CodeInterface codeInterface = null;
        PoolInfo poolInfo = null;

        try {
            logger.info("Starting worker " + id + " running " + codeName);

            poolInfo = new PoolInfo(id, nrOfWorkers);

            if (poolInfo.getRank() == 0) {
                codeInterface = new CodeInterface(id, codeName, codeDir,
                        amuseHome, jni, poolInfo.getHostnames());

                new Shutdown(codeInterface);

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
        }
    }

}
