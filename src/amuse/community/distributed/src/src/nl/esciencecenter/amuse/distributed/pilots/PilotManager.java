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

import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
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
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.filesystems.FileSystem;
import nl.esciencecenter.xenon.filesystems.Path;
import nl.esciencecenter.xenon.schedulers.JobDescription;
import nl.esciencecenter.xenon.schedulers.JobStatus;
import nl.esciencecenter.xenon.schedulers.Scheduler;
import nl.esciencecenter.xenon.utils.JavaJobDescription;

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

    private static int MAX_IO_LINES = 10000;

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
            int timeMinutes, int slotsPerNode, String nodeLabel, String options, String serverAddress, String[] hubAddresses,
            Path stdoutPath, Path stderrPath, boolean debug) throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        if (stdoutPath != null) {
            result.setStdout(stdoutPath.toRelativePath().toAbsolutePath().toString());
        }

        if (stderrPath != null) {
            result.setStderr(stderrPath.toRelativePath().toAbsolutePath().toString());
        }

        //~ result.setInteractive(false);

        if (queueName != null && !queueName.isEmpty()) {
            result.setQueueName(queueName);
        } else if (resource.isLocal() || resource.getSchedulerType().equals("ssh")) {
            result.setQueueName("unlimited");
        }

        if (resource.getSchedulerType().equals("slurm")) {
            result.setStartSingleProcess(true);
//            result.addJobOption("single.process", "true");
            result.setProcessesPerNode(slotsPerNode);

            //disable processor affinity
            result.addEnvironment("SLURM_CPU_BIND", "no");
            result.addEnvironment("SBATCH_CPU_BIND", "no");
        }

        //parse and add job options
        for (Map.Entry<String, String> option : parseOptions(options).entrySet()) {
            result.addJobOption(option.getKey(), option.getValue());
        }

        result.setNodeCount(nodeCount);
        result.setMaxRuntime(timeMinutes);

        AmuseConfiguration configuration = resource.getConfiguration();

        result.setExecutable(configuration.getJava());

        List<String> classpath = result.getJavaClasspath();
        classpath.add(configuration.getAmuseHome().getPath() + "/community/distributed/data/");
        classpath.add(configuration.getAmuseHome().getPath() + "/community/distributed/data/*");

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
    private final int slotsPerNode;
    private final int timeMinutes;
    private final int slots;
    private final String label;
    private final String options;

    private final String xenonJob;

    private final ResourceManager resource;

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
    public PilotManager(ResourceManager resource, String queueName, int nodeCount, int timeMinutes, int slotsPerNode, String nodeLabel,
            String options, String serverAddress, String[] hubAddresses, UUID amuseID, boolean debug)
            throws DistributedAmuseException {
        this.jobs = new ArrayList<AmuseJob>();
        ibisIdentifier = null;

        this.resource = resource;

        this.queueName = queueName;
        this.nodeCount = nodeCount;
        this.timeMinutes = timeMinutes;
        this.slotsPerNode = slotsPerNode;
        this.label = nodeLabel;
        this.options = options;

        this.slots = nodeCount * slotsPerNode;
        
        this.id = getNextID();

        try {
            Scheduler scheduler = resource.getScheduler();
            FileSystem filesystem = resource.getFileSystem();

            if (debug) {

                Path resourceHome = resource.getHome();

                //~ Path logDir = Utils.resolveWithRoot(xenon.files(), resourceHome, "distributed-amuse-logs");
                Path logDir = resourceHome.resolve("distributed-amuse-logs");

                logger.debug("logs will be put in dir: " + logDir);

                if (!filesystem.exists(logDir)) {
                    filesystem.createDirectories(logDir);
                }

                //~ stdoutPath = Utils.resolveWithRoot(xenon.files(), logDir, "amuse-" + amuseID + "-pilot-" + id + "-stdout.txt");
                //~ stderrPath = Utils.resolveWithRoot(xenon.files(), logDir, "amuse-" + amuseID + "-pilot-" + id + "-stderr.txt");
                stdoutPath = logDir.resolve("amuse-" + amuseID + "-pilot-" + id + "-stdout.txt");
                stderrPath = logDir.resolve("amuse-" + amuseID + "-pilot-" + id + "-stderr.txt");

            } else {
                stdoutPath = null;
                stderrPath = null;
            }

            //            Path xenonTmpDir = Utils.fromLocalPath(xenon.files(), tmpDir.getAbsolutePath());
            //            stdoutPath = Utils.resolveWithRoot(xenon.files(), xenonTmpDir, "reservation-" + uniqueID + "-stdout.txt");
            //            stderrPath = Utils.resolveWithRoot(xenon.files(), xenonTmpDir, "reservation-" + uniqueID + "-stderr.txt");

            JobDescription jobDescription = createJobDesciption(id, resource, queueName, nodeCount, timeMinutes, slotsPerNode,
                    nodeLabel, options, serverAddress, hubAddresses, stdoutPath, stderrPath, debug);

            logger.debug("starting reservation using scheduler {}", scheduler);

            this.xenonJob = scheduler.submitBatchJob(jobDescription);
            
            logger.debug("submitted reservation: {}", xenonJob);

            //get initial job status
            this.xenonJobStatus = scheduler.getJobStatus(this.xenonJob);
            
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

    public ResourceManager getResource() {
        return resource;
    }

    public int getResourceID() {
        return resource.getId();
    }

    public String getXenonJob() {
        return xenonJob;
    }
    
    private List<String> readlines(Path path) throws DistributedAmuseException{
        try {
            ArrayList<String> result = new ArrayList<>();
            InputStream stream = resource.getFileSystem().readFromFile(stdoutPath);
            BufferedReader buf= new BufferedReader( new InputStreamReader( stream, StandardCharsets.UTF_8) );
            String line = buf.readLine();
            while( line != null) {
              result.add( line );
              if(MAX_IO_LINES>0 && result.size() > MAX_IO_LINES) result.remove(0);
              line = buf.readLine();
            }
            stream.close();
            return result;
      } catch (Exception e) {
            throw new DistributedAmuseException("failed to read file " + path + " for " + this, e);
      }
    }
    
    public List<String> getStdout() throws DistributedAmuseException {
        if (stdoutPath == null) {
            ArrayList<String> result = new ArrayList<>();
            result.add("Disabled");
            return result;
        }
        try {
            return readlines(stdoutPath);            
        } catch (DistributedAmuseException e) {
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
            return readlines(stderrPath);            
        } catch (DistributedAmuseException e) {
            throw new DistributedAmuseException("failed to read stderr file " + stderrPath + " for " + this, e);
        }
    }

    public void stop() throws DistributedAmuseException {
        if (isDone()) {
            return;
        }

        logger.debug("cancelling xenon job for pilot: {}", this);
        try {
            resource.getScheduler().cancelJob(xenonJob);
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
        result.put("Slots Per Node", Integer.toString(slotsPerNode));
        result.put("Total Slots", Integer.toString(slots));
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
        if (!this.left && this.ibisIdentifier == null) {
            return "WAITING";
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
