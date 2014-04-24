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
package nl.esciencecenter.amuse.distributed.pilots;

import ibis.ipl.IbisIdentifier;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJob;
import nl.esciencecenter.amuse.distributed.remote.Pilot;
import nl.esciencecenter.amuse.distributed.resources.ResourceManager;
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.files.Path;
import nl.esciencecenter.xenon.jobs.Job;
import nl.esciencecenter.xenon.jobs.JobDescription;
import nl.esciencecenter.xenon.jobs.JobStatus;
import nl.esciencecenter.xenon.jobs.Scheduler;
import nl.esciencecenter.xenon.util.JavaJobDescription;
import nl.esciencecenter.xenon.util.Utils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Submits, monitors and controls pilot on some remote resource
 * 
 * 
 * @author Niels Drost
 * 
 */
public class PilotManager {

    public static final String WHITESPACE_REGEX = ";";

    public static final String EQUALS_REGEX = "\\s*=\\s*";

    private static final Logger logger = LoggerFactory.getLogger(PilotManager.class);

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

    private static JavaJobDescription createJobDesciption(int id, ResourceManager resource, String queueName, int nodeCount,
            int timeMinutes, int slots, String nodeLabel, String options, String serverAddress, String[] hubAddresses,
            Path stdoutPath, Path stderrPath, boolean debug) throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        if (stdoutPath != null) {
            result.setStdout(stdoutPath.getRelativePath().getAbsolutePath());
        }

        if (stderrPath != null) {
            result.setStderr(stderrPath.getRelativePath().getAbsolutePath());
        }

        result.setInteractive(false);

        if (queueName != null && !queueName.isEmpty()) {
            result.setQueueName(queueName);
        } else if (resource.isLocal() || resource.getSchedulerType().equals("ssh")) {
            result.setQueueName("unlimited");
        }

        if (resource.getSchedulerType().equals("slurm")) {
            result.addJobOption("single.process", "true");
            result.setProcessesPerNode(slots);

            //disable processor affinity
            result.addEnvironment("SLURM_CPU_BIND", "no");
            result.addEnvironment("SBATCH_CPU_BIND", "no");
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

        javaArguments.add("--pilot-id");
        javaArguments.add(Integer.toString(id));

        javaArguments.add("--resource-name");
        javaArguments.add(resource.getName());

        javaArguments.add("--server-address");
        javaArguments.add(serverAddress);

        javaArguments.add("--amuse-home");
        javaArguments.add(configuration.getAmuseHome().getAbsolutePath());
        
        if (resource.getTmpDir() != null) {
            javaArguments.add("--tmp-dir");
            javaArguments.add(resource.getTmpDir());
        }

        if (debug) {
            javaArguments.add("--debug");
        }

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

    private final String queueName;
    private final int nodeCount;
    private final int timeMinutes;
    private final int slots;
    private final String label;
    private final String options;

    private final Job xenonJob;

    private final ResourceManager resource;

    private final Xenon xenon;

    private final Path stdoutPath;
    private final Path stderrPath;

    private final List<AmuseJob> jobs;

    private IbisIdentifier ibisIdentifier;

    private JobStatus xenonJobStatus;

    private boolean left = false;
    
    private String lastState = null;

    /**
     * @param resource
     * @param queueName
     * @param nodeCount
     * @param timeMinutes
     * @param nodeLabel
     */
    public PilotManager(ResourceManager resource, String queueName, int nodeCount, int timeMinutes, int slots, String nodeLabel,
            String options, String serverAddress, String[] hubAddresses, Xenon xenon, File tmpDir, UUID amuseID, boolean debug)
            throws DistributedAmuseException {
        this.jobs = new ArrayList<AmuseJob>();
        ibisIdentifier = null;

        this.xenon = xenon;
        this.resource = resource;

        this.queueName = queueName;
        this.nodeCount = nodeCount;
        this.timeMinutes = timeMinutes;
        this.slots = slots;
        this.label = nodeLabel;
        this.options = options;

        this.id = getNextID();

        try {
            Scheduler scheduler = resource.getScheduler();

            if (debug) {

                Path resourceHome = resource.getHome();

                Path logDir = Utils.resolveWithRoot(xenon.files(), resourceHome, "distributed-amuse-logs");

                logger.debug("logs will be put in dir: " + logDir);

                if (!xenon.files().exists(logDir)) {
                    xenon.files().createDirectories(logDir);
                }

                stdoutPath = Utils.resolveWithRoot(xenon.files(), logDir, "amuse-" + amuseID + "-pilot-" + id + "-stdout.txt");
                stderrPath = Utils.resolveWithRoot(xenon.files(), logDir, "amuse-" + amuseID + "-pilot-" + id + "-stderr.txt");

            } else {
                stdoutPath = null;
                stderrPath = null;
            }

            //            Path xenonTmpDir = Utils.fromLocalPath(xenon.files(), tmpDir.getAbsolutePath());
            //            stdoutPath = Utils.resolveWithRoot(xenon.files(), xenonTmpDir, "reservation-" + uniqueID + "-stdout.txt");
            //            stderrPath = Utils.resolveWithRoot(xenon.files(), xenonTmpDir, "reservation-" + uniqueID + "-stderr.txt");

            JobDescription jobDescription = createJobDesciption(id, resource, queueName, nodeCount, timeMinutes, slots,
                    nodeLabel, options, serverAddress, hubAddresses, stdoutPath, stderrPath, debug);

            logger.debug("starting reservation using scheduler {}", scheduler);

            this.xenonJob = xenon.jobs().submitJob(scheduler, jobDescription);
            
            logger.debug("submitted reservation: {}", xenonJob);

            //get initial job status
            this.xenonJobStatus = xenon.jobs().getJobStatus(this.xenonJob);
            
            logState();

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

    public String getLabel() {
        return label;
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

    public Job getXenonJob() {
        return xenonJob;
    }
    
    public List<String> getStdout() throws DistributedAmuseException {
        if (stdoutPath == null) {
            ArrayList<String> result = new ArrayList<>();
            result.add("Disabled");
            return result;
        }
        try {
            return Utils.readAllLines(xenon.files(), stdoutPath, StandardCharsets.UTF_8);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to read stdout file " + stdoutPath + " for " + this, e);
        }
    }

    public List<String> getStderr() throws DistributedAmuseException {
        if (stderrPath == null) {
            ArrayList<String> result = new ArrayList<>();
            result.add("Disabled");
            return result;
        }
        try {
            return Utils.readAllLines(xenon.files(), stderrPath, StandardCharsets.UTF_8);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to read stderr file " + stderrPath + " for " + this, e);
        }
    }

    public void stop() throws DistributedAmuseException {
        if (isDone()) {
            return;
        }

        logger.debug("cancelling xenon job for pilot: {}", this);
        try {
            xenon.jobs().cancelJob(xenonJob);
        } catch (XenonException e) {
            throw new DistributedAmuseException("failed to cancel job " + xenonJob, e);
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
        return "Pilot [id=" + id + "]";
    }

    public Map<String, String> getStatusMap() {
        Map<String, String> result = new LinkedHashMap<String, String>();

        result.put("ID", Integer.toString(id));
        result.put("Queue", queueName);
        result.put("Node Count", Integer.toString(nodeCount));
        result.put("Time(minutes)", Integer.toString(timeMinutes));
        result.put("Slots", Integer.toString(slots));
        result.put("Node Label", label);
        result.put("Resource Name", getResourceName());
        result.put("Resource ID", Integer.toString(getResourceID()));

        result.put("Ibis Identifier", String.valueOf(getIbisIdentifier()));
        result.put("Running", Boolean.toString(isRunning()));
        result.put("Left", Boolean.toString(hasLeft()));
        result.put("Done", Boolean.toString(isDone()));

        result.put("Xenon Job Status", String.valueOf(getXenonJobStatus()));

        return result;
    }
    
    
    private synchronized void logState() {
        String state = getStateString();
        
        if (!state.equals(this.lastState)) {
            this.lastState = state;
            logger.info("State for pilot {} on resource {} now {}", getID(), getResourceName(), getStateString());    
        }
    }
    
    //ibis identifier, set by status monitor
    synchronized void setIbisIdentifier(IbisIdentifier ibis) {
        this.ibisIdentifier = ibis;
        notifyAll();
        logState();
    }

    public synchronized IbisIdentifier getIbisIdentifier() {
        return this.ibisIdentifier;
    }

    //status, set by xenon job monitor
    synchronized void setXenonJobStatus(JobStatus status) {
        this.xenonJobStatus = status;
        notifyAll();
        logState();
    }

    //status, set by Xenon job monitor
    synchronized JobStatus getXenonJobStatus() {
        return this.xenonJobStatus;
    }

    //the pilot has left
    synchronized void setLeft() {
        this.left = true;
        notifyAll();
        logState();
    }

    synchronized boolean hasLeft() {
        return this.left;
    }

    public synchronized boolean isRunning() {
        return !this.left && this.ibisIdentifier != null;
    }

    public synchronized boolean hasException() {
        return xenonJobStatus.hasException();
    }

    public synchronized boolean isDone() {
        //done if xenon job is done and either this pilot has left, or it has never started in the first place.
        return xenonJobStatus.isDone() && (this.left || this.ibisIdentifier == null) ;
    }

    public synchronized Exception getException() {
        return xenonJobStatus.getException();
    }

    public synchronized String getStateString() {
        if (isRunning()) {
            return "RUNNING";
        }

        return xenonJobStatus.getState();
    }

    public synchronized void addAmuseJob(AmuseJob amuseJob) {
        jobs.add(amuseJob);
    }

    private synchronized int availableSlots() {
        int usedSlotCount = 0;

        Iterator<AmuseJob> iterator = jobs.iterator();

        //ignore jobs that have already finished
        while (iterator.hasNext()) {
            AmuseJob job = iterator.next();

            if (job.isDone()) {
                iterator.remove();
            } else {
                usedSlotCount += job.getDescription().getNrOfSlots();
            }
        }

        int result = slots - usedSlotCount;

        logger.debug("pilot {} currently running {} jobs, has {} slots available", this, jobs.size(), result);

        return result;
    }

    public boolean canRun(AmuseJob job) {
        //check if we are running
        if (!isRunning()) {
            return false;
        }

        //check if the label matches (if given)
        if (job.getDescription().getLabel() != null && !job.getDescription().getLabel().equals(this.label)) {
            return false;
        }

        //check if we have enough slots available
        return job.getDescription().getNrOfSlots() <= availableSlots();
    }

}
