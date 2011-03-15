package ibis.amuse;

import ibis.ipl.Ibis;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReceivePort;
import ibis.ipl.RegistryEventHandler;
import ibis.ipl.SendPort;

import java.io.IOException;
import java.util.UUID;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class RemoteWorker implements RegistryEventHandler {
    
    private static class Shutdown extends Thread {
        private final RemoteWorker worker;

        Shutdown(RemoteWorker worker) {
            this.worker = worker;
        }

        public void run() {
            worker.end();
        }
    }
    
    private static Logger logger = LoggerFactory.getLogger(RemoteWorker.class);

    private final CommunityCode communityCode;

    private final SendPort sendPort;

    private final ReceivePort receivePort;

    private final Ibis ibis;

    private volatile boolean ended;

    public RemoteWorker(UUID id, String codeName) throws Exception {
        ibis = IbisFactory.createIbis(Daemon.ibisCapabilities, this,
                Daemon.portType);

        sendPort = ibis.createSendPort(Daemon.portType);
        receivePort = ibis.createReceivePort(Daemon.portType, id.toString());

        communityCode = new CommunityCode(codeName, receivePort, sendPort);

        receivePort.enableConnections();
        ibis.registry().enableEvents();
    }

    public static void main(String[] arguments) {

        String idString = System.getProperty("ibis.deploy.job.id");
        String codeName = null;

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equals("-i") || arguments[i].equals("--worker-id")) {
                i++;
                idString = arguments[i];
            } else if (arguments[i].equals("-c")
                    || arguments[i].equals("--code-name")) {
                i++;
                codeName = arguments[i];
            }
        }

        UUID id = UUID.fromString(idString);

        try {

            RemoteWorker worker = new RemoteWorker(id, codeName);
            
            new Shutdown(worker);

            worker.run();

        } catch (Exception e) {
            logger.error("Error on running remote worker", e);
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
            receivePort.close();
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
        end();
    }

}
