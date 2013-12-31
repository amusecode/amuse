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
package nl.esciencecenter.amuse.distributed.reservations;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.pilot.Pilot;
import nl.esciencecenter.amuse.distributed.resources.Resource;
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.adaptors.ssh.SshAdaptor;
import nl.esciencecenter.xenon.credentials.Credential;
import nl.esciencecenter.xenon.files.Path;
import nl.esciencecenter.xenon.jobs.Job;
import nl.esciencecenter.xenon.jobs.JobDescription;
import nl.esciencecenter.xenon.jobs.Scheduler;
import nl.esciencecenter.xenon.util.JavaJobDescription;
import nl.esciencecenter.xenon.util.Utils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class Reservation {

    public static final String WHITESPACE_REGEX = ";";

    public static final String EQUALS_REGEX = "\\s*=\\s*";

    private static final Logger logger = LoggerFactory.getLogger(Reservation.class);

    private static int nextID = 0;

    private static int getNextID() {
        return nextID++;
    }

    private static Map<String, String> parseOptions(String options) throws DistributedAmuseException {
        Map<String, String> result = new HashMap<String, String>();

        if (options == null || options.isEmpty()) {
            return result;
        }

        for (String option : options.split(WHITESPACE_REGEX)) {
            String[] keyvalue = option.split(EQUALS_REGEX, 2);
            if (keyvalue.length != 2) {
                throw new DistributedAmuseException("Key-Value option " + "\"" + option + "\" not a valid key=value pair");
            }
            logger.debug("adding option \"{}\" = \"{}\"", keyvalue[0], keyvalue[1]);
            result.put(keyvalue[0], keyvalue[1]);
        }

        return result;
    }

    private static JavaJobDescription createJobDesciption(int id, UUID uniqueID, Resource resource, String queueName,
            int nodeCount, int timeMinutes, int slots, String nodeLabel, String options, String serverAddress,
            String[] hubAddresses, Path stdoutPath, Path stderrPath) throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        result.setStdout(stdoutPath.getRelativePath().getAbsolutePath());
        result.setStderr(stderrPath.getRelativePath().getAbsolutePath());

        result.setInteractive(false);

        if (queueName != null && !queueName.isEmpty()) {
            result.setQueueName(queueName);
        }

        //parse and add job options
        for (Map.Entry<String, String> option : parseOptions(options).entrySet()) {
            result.addJobOption(option.getKey(), option.getValue());
        }

        result.setNodeCount(nodeCount);
        result.setMaxTime(timeMinutes);

        AmuseConfiguration configuration = resource.getConfiguration();

        result.setExecutable(configuration.getJava());

        List<String> classpath = result.getJavaClasspath();
        classpath.add(configuration.getAmuseHome().getPath() + "/src/amuse/community/distributed/src/dist/*");
        classpath.add(configuration.getAmuseHome().getPath() + "/src/amuse/community/distributed/worker.jar");
        classpath.add(configuration.getAmuseHome().getPath() + "/src/amuse/community/distributed");

        result.setJavaMain(Pilot.class.getCanonicalName());

        List<String> javaArguments = result.getJavaArguments();

        javaArguments.add("--reservation-id");
        javaArguments.add(Integer.toString(id));

        javaArguments.add("--node-label");
        javaArguments.add(nodeLabel);

        javaArguments.add("--resource-name");
        javaArguments.add(resource.getName());

        javaArguments.add("--server-address");
        javaArguments.add(serverAddress);

        javaArguments.add("--amuse-home");
        javaArguments.add(configuration.getAmuseHome().getAbsolutePath());

        javaArguments.add("--slots");
        javaArguments.add(Integer.toString(slots));

        String hubs = null;

        if (resource.getHub() != null) {
            hubs = resource.getHub().getAddress();
        } else {
            for (String hub : hubAddresses) {
                if (hubs == null) {
                    hubs = hub;
                } else {
                    hubs = hubs + "," + hub;
                }
            }
        }

        if (hubs != null) {
            javaArguments.add("--hub-addresses");
            javaArguments.add(hubs);
        }

        return result;
    }

    //ID in amuse
    private final int id;

    //unique ID for log files and such
    private final UUID uniqueID;

    private final String queueName;
    private final int nodeCount;
    private final int timeMinutes;
    private final int slots;
    private final String nodeLabel;
    private final String options;

    private final Job job;

    private final Resource resource;

    private final Xenon xenon;

    private final Path stdoutPath;
    private final Path stderrPath;

    /**
     * @param resource
     * @param queueName
     * @param nodeCount
     * @param timeMinutes
     * @param nodeLabel
     */
    public Reservation(Resource resource, String queueName, int nodeCount, int timeMinutes, int slots, String nodeLabel,
            String options, String serverAddress, String[] hubAddresses, Xenon xenon, File tmpDir)
            throws DistributedAmuseException {
        this.xenon = xenon;
        this.resource = resource;

        this.queueName = queueName;
        this.nodeCount = nodeCount;
        this.timeMinutes = timeMinutes;
        this.slots = slots;
        this.nodeLabel = nodeLabel;
        this.options = options;

        this.id = getNextID();
        this.uniqueID = UUID.randomUUID();

        try {
            Scheduler scheduler = resource.getScheduler();
            Path resourceHome = resource.getHome();

            Path logDir = Utils.resolveWithRoot(xenon.files(), resourceHome, "distributed-amuse-logs");
            
            logger.debug("logs will be put in dir: " + logDir);

            if (!xenon.files().exists(logDir)) { 
                xenon.files().createDirectories(logDir);
            }

            stdoutPath = Utils.resolveWithRoot(xenon.files(), logDir, "reservation-" + uniqueID + "-stdout.txt");
            stderrPath = Utils.resolveWithRoot(xenon.files(), logDir, "reservation-" + uniqueID + "-stderr.txt");
 
//            Path xenonTmpDir = Utils.fromLocalPath(xenon.files(), tmpDir.getAbsolutePath());
//            stdoutPath = Utils.resolveWithRoot(xenon.files(), xenonTmpDir, "reservation-" + uniqueID + "-stdout.txt");
//            stderrPath = Utils.resolveWithRoot(xenon.files(), xenonTmpDir, "reservation-" + uniqueID + "-stderr.txt");

            JobDescription jobDescription = createJobDesciption(id, uniqueID, resource, queueName, nodeCount, timeMinutes, slots,
                    nodeLabel, options, serverAddress, hubAddresses, stdoutPath, stderrPath);

            logger.debug("starting reservation using scheduler {}", scheduler);

            this.job = xenon.jobs().submitJob(scheduler, jobDescription);

            logger.debug("submitted reservation: {}", job);

        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start reservation on " + resource.getName() + ": " + e, e);
        }
    }

    public int getID() {
        return id;
    }

    public String getQueueName() {
        return queueName;
    }

    public int getNodeCount() {
        return nodeCount;
    }

    public int getTimeMinutes() {
        return timeMinutes;
    }

    public int getSlotsPerNode() {
        return slots;
    }

    public String getNodeLabel() {
        return nodeLabel;
    }
    
    public String getOptions() {
        return options;
    }

    public String getResourceName() {
        return resource.getName();
    }

    public int getResourceID() {
        return resource.getId();
    }

    public Job getJob() {
        return job;
    }

    public List<String> getStdout() throws DistributedAmuseException {
        try {
            return Utils.readAllLines(xenon.files(), stdoutPath, StandardCharsets.UTF_8);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to read stdout file " + stdoutPath + " for " + this, e);
        }
    }

    public List<String> getStderr() throws DistributedAmuseException {
        try {
            return Utils.readAllLines(xenon.files(), stderrPath, StandardCharsets.UTF_8);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to read stderr file " + stderrPath + " for " + this, e);
        }
    }

    public void cancel() throws DistributedAmuseException {
        logger.debug("cancelling reservation: {}", this);
        try {
            xenon.jobs().cancelJob(job);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to cancel job " + job, e);
        }
    }

    //    public void waitUntilStarted() throws DistributedAmuseException {
    //        try {
    //            JobStatus status = Xenon.jobs().waitUntilRunning(job, 0);
    //
    //            if (status.hasException()) {
    //                throw new DistributedAmuseException("error in reservation job: " + job, status.getException());
    //            }
    //        } catch (XenonIOException | XenonException e) {
    //            throw new DistributedAmuseException("failed to get job status " + job, e);
    //        }
    //    }
    //
    //    public String getStatus() throws DistributedAmuseException {
    //        try {
    //            return Xenon.jobs().getJobStatus(job).getState();
    //        } catch (XenonIOException | XenonException e) {
    //            throw new DistributedAmuseException("failed to get job status " + job, e);
    //        }
    //    }

    @Override
    public String toString() {
        return "Reservation [id=" + id + "]";
    }

    public Map<String, String> getStatusMap() {
        Map<String, String> result = new LinkedHashMap<String, String>();

        result.put("ID", Integer.toString(id));
        result.put("Queue", queueName);
        result.put("Node Count", Integer.toString(nodeCount));
        result.put("Time(minutes)", Integer.toString(timeMinutes));
        result.put("Slots", Integer.toString(slots));
        result.put("Node Label", nodeLabel);
        result.put("Resource Name", getResourceName());
        result.put("Resource ID", Integer.toString(getResourceID()));

        return result;
    }

}
