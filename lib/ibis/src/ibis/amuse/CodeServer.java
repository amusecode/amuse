package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReceivePort;
import ibis.ipl.RegistryEventHandler;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CodeServer implements RegistryEventHandler {

    private static class Shutdown extends Thread {
        private final CodeServer worker;

        Shutdown(CodeServer worker) {
            this.worker = worker;
        }

        public void run() {
            worker.end();
        }
    }

    private static Logger logger = LoggerFactory.getLogger(CodeServer.class);

    private final Runnable communityCode;

    private final SendPort sendPort;

    private final ReceivePort receivePort;

    private final Ibis ibis;

    private volatile boolean ended = false;

    public CodeServer(String id, String codeName, String codeDir, String amuseHome, boolean jni) throws Exception {
        ibis = IbisFactory.createIbis(Daemon.ibisCapabilities, this,
                Daemon.portType);

        sendPort = ibis.createSendPort(Daemon.portType);
        receivePort = ibis.createReceivePort(Daemon.portType, id);

        if (jni) {
            communityCode = new JNICode(codeName, receivePort, sendPort);
        } else {
            communityCode = new SocketCode(codeName, codeDir, receivePort, sendPort);  
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

    public void run() {

        communityCode.run();

        end();
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
    
    public static void main(String[] arguments) {

        String id = null;
        String codeName = null;
        boolean jni = false;
        String amuseHome = null;
        String codeDir = null;

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
            }
        }

        try {
            logger.info("Starting worker " + id + " running " + codeName);

            CodeServer server = new CodeServer(id, codeName, codeDir, amuseHome, jni);

            new Shutdown(server);

            server.run();
            logger.info("Worker " + id + " running " + codeName + " DONE!");
        } catch (Exception e) {
            logger.error("Error on running remote worker", e);
        }
    }


}
