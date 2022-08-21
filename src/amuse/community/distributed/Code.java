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
import java.util.List;
import java.util.UUID;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJob;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJobDescription;
import nl.esciencecenter.amuse.distributed.jobs.FunctionJob;
import nl.esciencecenter.amuse.distributed.jobs.FunctionJobDescription;
import nl.esciencecenter.amuse.distributed.jobs.ScriptJob;
import nl.esciencecenter.amuse.distributed.jobs.ScriptJobDescription;
import nl.esciencecenter.amuse.distributed.jobs.WorkerJobDescription;
import nl.esciencecenter.amuse.distributed.jobs.WorkerJob;
import nl.esciencecenter.amuse.distributed.jobs.WorkerJob;
import nl.esciencecenter.amuse.distributed.pilots.PilotManager;
import nl.esciencecenter.amuse.distributed.resources.ResourceManager;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Interface from generated code interface to the distributed amuse implementation. Simply forwards all calls. Must be in the
 * "default" package as the CodeInterface and Worker are also generated there.
 * 
 * @author Niels Drost
 * 
 */
public class Code implements CodeInterface {

    private static final Logger logger = LoggerFactory.getLogger(Code.class);

    private final String codeDir;

    private final String amuseRootDir;

    private DistributedAmuse distributedAmuse;

    private String currentError = "";

    private boolean debug = false;
    private boolean startHubs = true;
    private int webinterfacePort = 0;
    private int workerStartupTimeout = 60;
    private int workerQueueTimeout = 60;

    private final UUID id;

    public Code(String codeDir, String amuseRootDir) throws DistributedAmuseException {
        this.codeDir = codeDir;
        this.amuseRootDir = amuseRootDir;

        id = UUID.randomUUID();

        distributedAmuse = null;
    }

    @Override
    public synchronized int get_current_error(String[] result) {
        result[0] = currentError;
        return 0;
    }

    private synchronized void reportError(String message, DistributedAmuseException e) {
        logger.error(message, e);
        currentError = message;
    }

    public int initialize_code() {
        return 0;
    }

    @Override
    public int commit_parameters() {
        if (distributedAmuse == null) {
            try {
                distributedAmuse = new DistributedAmuse(id, codeDir, amuseRootDir, webinterfacePort, debug, startHubs,
                        workerQueueTimeout, workerStartupTimeout);
            } catch (DistributedAmuseException e) {
                logger.error("Exception while initializing code", e);
                return -10;
            }
            return 0;
        } else {
            return -10;
        }
    }

    @Override
    public int recommit_parameters() {
        logger.error("Recommit of parameters not (yet) supported");
        return -1;
    }

    @Override
    public int get_debug(int[] debug) {
        debug[0] = booleanToInteger(this.debug);
        return 0;
    }

    @Override
    public int set_debug(int debug) {
        this.debug = integerToBoolean(debug);
        return 0;
    }

    @Override
    public int get_worker_port(int[] result) {
        logger.debug("Returning worker port.");
        result[0] = -1;
        if(distributedAmuse != null)  result[0] = distributedAmuse.getWorkerPort();

        return 0;
    }

    @Override
    public int get_start_hubs(int[] start_hubs) {
        start_hubs[0] = booleanToInteger(this.startHubs);
        return 0;
    }

    @Override
    public int set_start_hubs(int start_hubs) {
        this.startHubs = integerToBoolean(start_hubs);
        return 0;
    }

    @Override
    public int get_worker_startup_timeout(int[] result) {
        logger.debug("Returning worker startup timeout.");

        result[0] = this.workerStartupTimeout;

        return 0;
    }

    @Override
    public int set_worker_startup_timeout(int worker_startup_timeout) {
        this.workerStartupTimeout = worker_startup_timeout;
        return 0;
    }
    
    @Override
    public int get_worker_queue_timeout(int[] result) {
        logger.debug("Returning worker port.");

        result[0] = this.workerQueueTimeout;

        return 0;
    }

    @Override
    public int set_worker_queue_timeout(int worker_queue_timeout) {
        this.workerQueueTimeout = worker_queue_timeout;
        return 0;
    }


    @Override
    public int get_webinterface_port(int[] result) {
        logger.debug("Returning worker port.");

        if (distributedAmuse == null) {
            //not initialized yet, return current value of parameter
            result[0] = this.webinterfacePort;
        } else {
            //return resulting port
            result[0] = distributedAmuse.webInterface().getPort();
        }

        return 0;
    }

    @Override
    public int set_webinterface_port(int webinterface_port) {
        this.webinterfacePort = webinterface_port;
        return 0;
    }

    @Override
    public int new_resource(int[] resource_id, String[] name, String[] location, String[] amuse_dir, String[] tmp_dir, String[] gateway,
            String[] scheduler_type, String[] hub_queue_name, int[] hub_time_minutes,  int count) {
        try {
            for (int i = 0; i < count; i++) {
                String locationString = location[i];
                
                if (locationString.isEmpty()) {
                    locationString = null;
                }
                
                String gatewayString = gateway[i];
                
                if (gatewayString.isEmpty()) {
                    gatewayString = null;
                }
                
                String tmpDirString = tmp_dir[i];
                
                if (tmpDirString.isEmpty()) {
                    tmpDirString = null;
                }
                
                String hubQueueNameString = hub_queue_name[i];
                
                if (hubQueueNameString.isEmpty()) {
                    hubQueueNameString = null;
                }
                
                ResourceManager resource = distributedAmuse.resources().newResource(name[i], locationString, gatewayString,
                        amuse_dir[i], tmpDirString, scheduler_type[i], hubQueueNameString, hub_time_minutes[i]);
                resource_id[i] = resource.getId();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on creating new resource: " + e, e);
            return -10;
        }
    }

    //boolean as an int
    private static int booleanToInteger(boolean value) {
        if (value) {
            return 1;
        } else {
            return 0;
        }
    }

    private static boolean integerToBoolean(int value) {
        if (value <= 0) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public int get_resource_state(int[] index_of_the_resource, String[] name, String[] location, String[] gateway,
            String[] amuse_dir, String[] tmp_dir, String[] scheduler_type, String[] hub_queue_name, int[] hub_time_minutes, int count) {
        try {
            for (int i = 0; i < count; i++) {
                ResourceManager resource = distributedAmuse.resources().getResource(index_of_the_resource[i]);

                name[i] = resource.getName();
                location[i] = resource.getLocation();
                gateway[i] = resource.getGateway();
                amuse_dir[i] = resource.getAmuseDir();
                tmp_dir[i] = resource.getTmpDir();
                scheduler_type[i] = resource.getSchedulerType();
                hub_queue_name[i] = resource.getHubQueueName();
                hub_time_minutes[i] = resource.getHubTimeMinutes();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int delete_resource(int[] index_of_the_resource, int count) {
        try {
            for (int i = 0; i < count; i++) {
                ResourceManager resource = distributedAmuse.resources().getResource(index_of_the_resource[i]);
                distributedAmuse.resources().deleteResource(resource);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on deleting resource: " + e, e);
            return -10;
        }
    }

    @Override
    public int new_pilot(int[] pilot_id, String[] resource_name, String[] queue_name, int[] node_count, int[] time_minutes,
            int[] slots, String[] label, String[] options, int count) {
        try {
            for (int i = 0; i < count; i++) {
                PilotManager result = distributedAmuse.pilots().newPilot(resource_name[i], queue_name[i], node_count[i],
                        time_minutes[i], slots[i], label[i], options[i]);

                pilot_id[i] = result.getID();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on creating new pilot: " + e, e);
            return -10;
        }

    }

    @Override
    public int get_pilot_state(int[] pilot_id, String[] resource_name, String[] queue_name, int[] node_count, int[] time,
            int[] slots_per_node, String[] node_label, String[] status, String[] options, int count) {
        try {
            for (int i = 0; i < count; i++) {
                PilotManager pilot = distributedAmuse.pilots().getPilot(pilot_id[i]);

                resource_name[i] = pilot.getResourceName();
                queue_name[i] = pilot.getQueueName();
                node_count[i] = pilot.getNodeCount();
                time[i] = pilot.getTimeMinutes();
                slots_per_node[i] = pilot.getSlotsPerNode();
                node_label[i] = pilot.getLabel();
                options[i] = pilot.getOptions();
                status[i] = pilot.getStateString();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int get_pilot_status(int[] pilot_id, String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                PilotManager pilot = distributedAmuse.pilots().getPilot(pilot_id[i]);

                status[i] = pilot.getStateString();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int delete_pilot(int[] pilot_id, int count) {
        try {
            for (int i = 0; i < count; i++) {
                distributedAmuse.pilots().deletePilot(pilot_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on deleting pilot: " + e, e);
            return -10;
        }
    }

    @Override
    public int wait_for_pilots() {
        try {
            distributedAmuse.pilots().waitForAllPilots();
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on waiting for pilots: " + e, e);
            return -10;
        }
    }

    @Override
    public int submit_script_job(int[] job_id, String[] script_dir, String[] script_name, String[] arguments, String[] input_dir, String[] output_dir, String[] stdout_file, String[] stderr_file, String[] label, int count) {
        try {
            for (int i = 0; i < count; i++) {
                String stdoutFile = stdout_file[i];

                if (stdoutFile.isEmpty()) {
                    stdoutFile = null;
                }

                String stderrFile = stderr_file[i];

                if (stderrFile.isEmpty()) {
                    stderrFile = null;
                }

                String theLabel = label[i];

                if (theLabel.equals("")) {
                    theLabel = null;
                }
                
                String inputDir = input_dir[i];
                if (inputDir.isEmpty()) {
                    inputDir = null;
                }
                
                String outputDir = output_dir[i];
                if (outputDir.isEmpty()) {
                    outputDir = null;
                }


                ScriptJobDescription description = new ScriptJobDescription(stdoutFile, stderrFile, theLabel, script_name[i],
                        arguments[i], script_dir[i], inputDir, outputDir);

                ScriptJob job = distributedAmuse.jobs().submitScriptJob(description);
                job_id[i] = job.getJobID();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Cannot submit script job: " + e, e);
            return -10;
        }

    }

    @Override
    public int get_script_job_state(int[] job_id, String[] script_dir, String[] script_name, String[] arguments, String[] input_dir, String[] output_dir, String[] stdout_file, String[] stderr_file, String[] label, String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                ScriptJob job = distributedAmuse.jobs().getScriptJob(job_id[i]);

                ScriptJobDescription description = job.getDescription();

                script_dir[i] = description.getScriptDir();
                script_name[i] = description.getScriptName();
                arguments[i] = description.getArguments();
                input_dir[i] = description.getInputDir();
                output_dir[i] = description.getOutputDir();
                stdout_file[i] = description.getStdoutFile();
                stderr_file[i] = description.getStderrFile();
                label[i] = description.getLabel();

                status[i] = job.getJobState();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int get_script_job_status(int[] job_id, String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                ScriptJob job = distributedAmuse.jobs().getScriptJob(job_id[i]);

                status[i] = job.getJobState();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int delete_script_job(int[] job_id, int count) {
        try {
            for (int i = 0; i < count; i++) {
                //Delete from manager list. Will also cancel job
                distributedAmuse.jobs().removeJob(job_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on deleting pilot: " + e, e);
            return -10;
        }
    }

    @Override
    public int wait_for_script_jobs() {
        try {
            distributedAmuse.jobs().waitForScriptJobs();
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on waiting for jobs: " + e, e);
            return -10;
        }

    }

    @Override
    public int submit_function_job(int[] job_id, String[] function, String[] arguments, String[] kwarguments, String[] stdout_file,
            String[] stderr_file, String[] label, int count) {
        try {
            for (int i = 0; i < count; i++) {
                FunctionJobDescription description = new FunctionJobDescription(function[i], arguments[i], kwarguments[i], stdout_file[i],
                        stderr_file[i], label[i]);

                FunctionJob job = distributedAmuse.jobs().submitFunctionJob(description);

                job_id[i] = job.getJobID();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Cannot submit pickled job: " + e, e);
            return -10;
        }

    }

    @Override
    public int get_function_job_status(int[] job_id, String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                FunctionJob job = distributedAmuse.jobs().getFunctionJob(job_id[i]);

                status[i] = job.getJobState();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int get_function_job_state(int[] job_id, String[] stdout_file, String[] stderr_file, String[] node_label,
            String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                FunctionJob job = distributedAmuse.jobs().getFunctionJob(job_id[i]);

                AmuseJobDescription description = job.getDescription();

                stdout_file[i] = description.getStdoutFile();
                stderr_file[i] = description.getStderrFile();
                node_label[i] = description.getLabel();
                status[i] = job.getJobState();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int get_function_job_result(int[] job_id, String[] result, int count) {
        try {
            for (int i = 0; i < count; i++) {
                FunctionJob job = distributedAmuse.jobs().getFunctionJob(job_id[i]);
                result[i] = job.getResult();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting job result: " + e, e);
            return -10;
        }
    }

    @Override
    public int delete_function_job(int[] job_id, int count) {
        try {
            for (int i = 0; i < count; i++) {
                //delete from manager list, will also cancel job if needed
                distributedAmuse.jobs().removeJob(job_id[i]);
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on deleting function job: " + e, e);
            return -10;
        }
    }

    public int get_number_of_workers(int[] number_of_workers) {
        number_of_workers[0] = distributedAmuse.jobs().getWorkerIDs().size();
        return 0;
    }

    @Override
    /**
     * @param index only here to force AMUSE to generate the count parameter
     * @param id_of_the_worker outgoing indexes
     * @param count number of indexes to return
     * @return
     */
    public int get_worker_ids(int[] index, int[] id_of_the_worker, int count) {
        List<Integer> workerIDs = distributedAmuse.jobs().getWorkerIDs();
        
        if (workerIDs.size() != count) {
            reportError("Error on getting worker IDs: number of indexes requested " + count
                    + " does not match number of workers " + workerIDs.size(), null);
            return -10;
        }
        
        for (int i = 0; i < count; i++) {
            id_of_the_worker[i] = workerIDs.get(i);
        }
        return 0;
    }

    @Override
    public int get_worker_status(int[] worker_id, String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                WorkerJob worker = distributedAmuse.jobs().getWorkerJob(worker_id[i]);

                status[i] = worker.getJobState();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int get_worker_state(int[] worker_id, String[] executable, String[] node_label, int[] worker_count,
            int[] thread_count, String[] status, int count) {
        try {
            for (int i = 0; i < count; i++) {
                WorkerJob worker = distributedAmuse.jobs().getWorkerJob(worker_id[i]);

                WorkerJobDescription description = worker.getDescription();

                executable[i] = description.getExecutable();
                node_label[i] = description.getLabel();
                worker_count[i] = description.getNrOfWorkers();
                thread_count[i] = description.getNrOfThreads();

                status[i] = worker.getJobState();
            }
            return 0;
        } catch (DistributedAmuseException e) {
            reportError("Error on getting resource state: " + e, e);
            return -10;
        }
    }

    @Override
    public int end_all() {
        if (distributedAmuse != null) {
            distributedAmuse.end();
        }
        return 0;
    }

    @Override
    public int startup_viz() {
        try {
          if (distributedAmuse != null) {
              distributedAmuse.startup_viz();
          }
          return 0;
        } catch (DistributedAmuseException e) {
          reportError("Error on starting up viz", e);
          return -1;
        }
    }

    @Override
    public void end() {
        if (distributedAmuse != null) {
            distributedAmuse.end();
        }
    }

}
