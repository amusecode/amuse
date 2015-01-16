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
package nl.esciencecenter.amuse.distributed.remote;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCreationFailedException;
import ibis.ipl.IbisFactory;
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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJobDescription;

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

    private final Watchdog watchdog;

    private final boolean debug;
    
    private final Path tmpDir;

    private static void initializeLogger(boolean debug) {
        if (debug) {
            ch.qos.logback.classic.Logger amuseLogger = (ch.qos.logback.classic.Logger) LoggerFactory
                    .getLogger("nl.esciencecenter.amuse");

            amuseLogger.setLevel(Level.DEBUG);

            logger.debug("DEBUG Enabled");
        }
    }
    
    Pilot(AmuseConfiguration configuration, Properties properties, int id, boolean debug, String tmpDir) throws IbisCreationFailedException,
            IOException, InterruptedException {
        this.configuration = configuration;
        this.debug = debug;
        jobs = new HashMap<Integer, JobRunner>();

        initializeLogger(debug);

        //ID of this pilot
        String tag = Integer.toString(id);

        logger.debug("Creating Ibis");

        watchdog = new Watchdog();

        ibis = IbisFactory.createIbis(DistributedAmuse.IPL_CAPABILITIES, properties, true, watchdog, null, tag,
                DistributedAmuse.ONE_TO_ONE_PORT_TYPE, DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

        logger.debug("Creating Receive port");

        receivePort = ibis.createReceivePort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE, PORT_NAME, this, this, null);
        
        Path tmpPath = Paths.get(tmpDir);
        
        Files.createDirectories(tmpPath);
        
        this.tmpDir = Files.createTempDirectory(tmpPath, "distributed-amuse-pilot-" + id + "-");
        
        logger.debug("Saving temporary files in " + this.tmpDir);

    }

    /**
     * Small run function. Waits until the "main" distributed amuse node declares it is time to go.
     * 
     * @throws IOException
     */
    private void run() throws IOException {
        receivePort.enableConnections();
        receivePort.enableMessageUpcalls();
        ibis.registry().enableEvents();

        logger.info("Pilot fully initialized, waiting for commands...");

        //Wait until the pool is terminated by the DistributedAmuse master node, or the master node leaves, or the
        //watchdog expires
        boolean cleanExit = watchdog.waitUntilExpired();

        logger.info("Pool terminated, ending pilot");

        removeFinishedJobs();
        for (JobRunner job : getJobRunners()) {
            logger.info("Ending job: " + job);
            job.interrupt();
            try {
                job.join(1000);
            } catch (InterruptedException e) {
                //IGNORE
            }
            logger.info("Job {} ended", job);
        }

        if (cleanExit) {
            ibis.end();
        }
    }

    private synchronized void addJobRunner(int jobID, JobRunner jobRunner) {
        jobs.put(jobID, jobRunner);
    }

    private synchronized JobRunner getJobRunner(int jobID) {
        return jobs.get(jobID);
    }

    private synchronized JobRunner[] getJobRunners() {
        return jobs.values().toArray(new JobRunner[0]);
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
            ReceivePortIdentifier resultPort = (ReceivePortIdentifier) readMessage.readObject();
            AmuseJobDescription description = (AmuseJobDescription) readMessage.readObject();
            
            logger.debug("Running job " + description);
            
            String type = description.getType();

            try {
                JobRunner jobRunner;

                switch (type) {
                case "worker":
                    jobRunner = new WorkerJobRunner(description, configuration, resultPort, ibis, tmpDir, debug, readMessage);

                    addJobRunner(description.getID(), jobRunner);
                    break;
                case "script":
                    jobRunner = new ScriptJobRunner(description, configuration, resultPort, ibis, tmpDir, debug, readMessage);

                    addJobRunner(description.getID(), jobRunner);
                    break;
                case "function":
                    jobRunner = new FunctionJobRunner(description, configuration, resultPort, ibis, tmpDir, debug, readMessage);

                    addJobRunner(description.getID(), jobRunner);
                    break;
                default:
                    throw new Exception("Unknown job type: " + type);

                }
            } catch (Exception e) {
                logger.error("Error starting job", e);
                replyException = new Exception("Error starting job:" + e, e);
            }
            logger.debug("Sending reply");

            //send reply
            SendPort sendPort = ibis.createSendPort(DistributedAmuse.ONE_TO_ONE_PORT_TYPE);

            sendPort.connect(replyPort, 60000, true);

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
                jobRunner.cancel();
            }
        } else if (command.equals("ping")) {
            readMessage.finish();
            watchdog.gotPing();
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

    public static void main(String[] arguments) throws Exception {
        File amuseHome = null;
        int pilotID = 0;
        boolean debug = false;
        String tmpDir = System.getProperty("java.io.tmpdir") + File.separator + "distributed-amuse-pilots-" + System.getProperty("user.name");
                

        Properties properties = new Properties();
        properties.put(IbisProperties.POOL_NAME, "amuse");
        properties.put(IbisProperties.SERVER_IS_HUB, "false");
        //properties.put("ibis.managementclient", "true");
        //properties.put("ibis.bytescount", "true");

        for (int i = 0; i < arguments.length; i++) {
            if (arguments[i].equalsIgnoreCase("--pilot-id")) {
                i++;
                pilotID = Integer.parseInt(arguments[i]);
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
            } else if (arguments[i].equalsIgnoreCase("--debug")) {
                debug = true;
            } else if (arguments[i].equalsIgnoreCase("--tmp-dir")) {
                i++;
                tmpDir = arguments[i];
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

        Pilot pilot = new Pilot(configuration, properties, pilotID, debug, tmpDir);

        pilot.run();

        logger.debug("Main pilot thread ended");
    }

}
