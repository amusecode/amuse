/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nl.esciencecenter.amuse.distributed.pilot;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCreationFailedException;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.IbisProperties;
import ibis.ipl.MessageUpcall;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.ReceivePortConnectUpcall;
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.SendPortIdentifier;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.net.InetAddress;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.jobs.WorkerDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.qos.logback.classic.Level;

/**
 * Pilot job. Started when a reservations is made.
 * 
 * @author Niels Drost
 * 
 */
public class Pilot implements MessageUpcall, ReceivePortConnectUpcall {
    
    public static final String WHITESPACE_REGEX = "\\s+";

    public static final String PORT_NAME = "pilot";

    private static final Logger logger = LoggerFactory.getLogger(Pilot.class);

    private final Ibis ibis;

    private final ReceivePort receivePort;

    private final HashMap<Integer, JobRunner> jobs;

    private final AmuseConfiguration configuration;
    
    private UUID id;
    
    private File tmpDir;

    private static File createTmpDir(UUID id) throws IOException {
        File systemTmpDir = new File(System.getProperty("java.io.tmpdir"));
        String userName = System.getProperty("user.name");

        if (!systemTmpDir.exists()) {
            throw new IOException("Java tmpdir does not exist " + systemTmpDir);
        }

        File result = new File(systemTmpDir, "distributed-amuse-" + userName + "/pilot-" + id.toString());
        result.mkdirs();

        return result;
    }
    
    private static void initializeLogger(boolean debug) {
        if (debug) {
            ch.qos.logback.classic.Logger amuseLogger = (ch.qos.logback.classic.Logger) LoggerFactory
                    .getLogger("nl.esciencecenter.amuse");

            amuseLogger.setLevel(Level.DEBUG);
            
            logger.debug("DEBUG Enabled");
        }
    }
    
    private static void runBootCommand(String command) throws IOException, InterruptedException {
        ProcessBuilder builder = new ProcessBuilder(command.split(WHITESPACE_REGEX));
        
        for (String key : builder.environment().keySet().toArray(new String[0])) {
            for (String blacklistedKey : WorkerProxy.ENVIRONMENT_BLACKLIST) {
                if (key.startsWith(blacklistedKey)) {
                    builder.environment().remove(key);
                    logger.info("removed " + key + " from environment");
                }
            }
        }
        
        builder.redirectError(Redirect.INHERIT);
        builder.redirectOutput(Redirect.INHERIT);
        
        logger.info("running boot command: {}", builder.command());
        
        Process bootProcess = builder.start();
        
        int exitcode = bootProcess.waitFor();
        
        if (exitcode != 0) {
            throw new IOException("Error (exit status " + exitcode + ") while running boot command: " + command);
        }
        
    }

    Pilot(AmuseConfiguration configuration, Properties properties, String reservationID, String nodeLabel, int slots, String bootCommand, boolean debug)
            throws IbisCreationFailedException, IOException, InterruptedException {
        this.configuration = configuration;
        jobs = new HashMap<Integer, JobRunner>();

        initializeLogger(debug);
        
        //reservation ID, label, slots, hostname
        String tag = reservationID + "," + nodeLabel + "," + slots + "," + InetAddress.getLocalHost().getHostAddress();

        ibis = IbisFactory.createIbis(DistributedAmuse.IPL_CAPABILITIES, properties, true, null, null, tag,
                DistributedAmuse.ONE_TO_ONE_PORT_TYPE, DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

        receivePort = ibis.createReceivePort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE, PORT_NAME, this, this, null);
        
        this.id = UUID.randomUUID();
        tmpDir = createTmpDir(id);
        
        if (bootCommand != null) {
            runBootCommand(bootCommand);
        }
        
    }

    /**
     * Small run function. Waits until the "main" distributed amuse node declares it is time to go.
     * 
     * @throws IOException
     */
    private void run() throws IOException {
        receivePort.enableConnections();
        receivePort.enableMessageUpcalls();
        
        //ibis.registry().enableEvents();

        //Wait until the pool is terminated by the DistributedAmuse master node
        //FIXME: no way to interrupt this wait.
        ibis.registry().waitUntilTerminated();

        logger.info("Pool terminated, ending pilot");

        ibis.end();
    }

    public static void main(String[] arguments) throws Exception {
        File amuseHome = null;
        String nodeLabel = "default";
        String reservationID = null;
        int slots = 1;
        String bootCommand = null;
        boolean debug = false;

        Properties properties = new Properties();
        properties.put(IbisProperties.POOL_NAME, "amuse");
        properties.put(IbisProperties.SERVER_IS_HUB, "false");
        //properties.put("ibis.managementclient", "true");
        //properties.put("ibis.bytescount", "true");

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equalsIgnoreCase("--reservation-id")) {
                i++;
                reservationID = arguments[i];
            } else if (arguments[i].equalsIgnoreCase("--node-label")) {
                i++;
                nodeLabel = arguments[i];
            } else if (arguments[i].equalsIgnoreCase("--resource-name")) {
                i++;
                properties.put(IbisProperties.LOCATION, arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--server-address")) {
                i++;
                properties.put(IbisProperties.SERVER_ADDRESS, arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--amuse-home")) {
                i++;
                amuseHome = new File(arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--hub-addresses")) {
                i++;
                properties.put(IbisProperties.HUB_ADDRESSES, arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--slots")) {
                i++;
                slots = Integer.parseInt(arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--boot-command")) {
                i++;
                bootCommand = arguments[i];
            } else if (arguments[i].equalsIgnoreCase("--debug")) {
                debug = true;
            } else {
                System.err.println("Unknown command line option: " + arguments[i]);
                System.exit(1);
            }
        }

        AmuseConfiguration configuration = new AmuseConfiguration(amuseHome);

        System.err.println("running Pilot using properties:");
        for (Entry<Object, Object> entry : properties.entrySet()) {
            System.err.println(entry.getKey() + " = " + entry.getValue());
        }
        
        Pilot pilot = new Pilot(configuration, properties, reservationID, nodeLabel, slots, bootCommand, debug);

        pilot.run();

        logger.debug("Main pilot thread ended");
    }

    private synchronized void addJobRunner(int jobID, JobRunner jobRunner) {
        jobs.put(jobID, jobRunner);
    }

    private synchronized JobRunner getJobRunner(int jobID) {
        return jobs.get(jobID);
    }

    private synchronized void removeFinishedJobs() {
        Iterator<Map.Entry<Integer, JobRunner>> iterator = jobs.entrySet().iterator();
        //remove all finished jobs
        while (iterator.hasNext()) {
            Map.Entry<Integer, JobRunner> entry = iterator.next();
            if (!entry.getValue().isAlive()) {
                iterator.remove();
            }
        }

    }

    /**
     * Handle incoming message from the amuse node
     */
    @Override
    public void upcall(ReadMessage readMessage) throws IOException, ClassNotFoundException {
        logger.debug("Received message upcall");

        Exception replyException = null;

        //run cleanup
        removeFinishedJobs();

        String command = readMessage.readString();

        logger.debug("Got command: " + command);

        if (command.equals("start")) {
            ReceivePortIdentifier replyPort = (ReceivePortIdentifier) readMessage.readObject();

            //details of job
            int jobID = readMessage.readInt();

            ReceivePortIdentifier resultPort = (ReceivePortIdentifier) readMessage.readObject();

            //hard coded worker job info
            WorkerDescription description = (WorkerDescription) readMessage.readObject();
            IbisIdentifier[] nodes = (IbisIdentifier[]) readMessage.readObject();

            readMessage.finish();

            //FIXME: transfer files etc

            try {
                JobRunner jobRunner = new JobRunner(jobID, description, configuration, nodes, resultPort, ibis, tmpDir);

                addJobRunner(jobID, jobRunner);

                //start a thread to run job.
                jobRunner.start();
            } catch (Exception e) {
                logger.error("Error starting job", e);
                replyException = new Exception("Error starting job:" + e, e);
            }
            logger.debug("Sending reply");

            //send reply
            SendPort sendPort = ibis.createSendPort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE);

            sendPort.connect(replyPort);

            WriteMessage reply = sendPort.newMessage();

            reply.writeObject(replyException);
            
            reply.finish();

            sendPort.close();

            logger.debug("Reply sent");
        } else if (command.equals("cancel")) {
            int jobID = readMessage.readInt();
            readMessage.finish();

            JobRunner jobRunner = getJobRunner(jobID);

            if (jobRunner != null) {
                //signal the thread it is time to cancel the job
                jobRunner.interrupt();
            }
        } else {
            logger.error("Failed to handle message, unknown command: " + command);
            
        }

    }

    @Override
    public boolean gotConnection(ReceivePort receiver, SendPortIdentifier applicant) {
        logger.debug("Got connection from {}", applicant);
        return true;
    }

    /**
     * @param receiver
     * @param origin
     * @param cause
     */
    @Override
    public void lostConnection(ReceivePort receiver, SendPortIdentifier origin, Throwable cause) {
        logger.debug("Lost connection to {}", origin);

    }

}
