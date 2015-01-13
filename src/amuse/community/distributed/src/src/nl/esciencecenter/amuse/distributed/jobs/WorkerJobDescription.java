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
package nl.esciencecenter.amuse.distributed.jobs;

import java.io.Serializable;

/**
 * Description of a worker.
 * 
 * 
 * @author Niels Drost
 */
public class WorkerJobDescription extends AmuseJobDescription implements Serializable {

    private static final long serialVersionUID = 1L;

    private final String executable;
    private final String workerDir;

    private final int nrOfWorkers;
    private final int nrOfThreads;
    private final boolean dynamicPythonCode;

    private final int startupTimeout;

    public WorkerJobDescription(String stdoutFile, String stderrFile, String nodeLabel, String executable, String workerDir,
            int nrOfWorkers, int nrOfThreads, boolean dynamicPythonCode, int startupTimeout) {
        super(stdoutFile, stderrFile, nodeLabel);
        this.executable = executable;
        this.workerDir = workerDir;
        this.nrOfWorkers = nrOfWorkers;
        this.nrOfThreads = nrOfThreads;
        this.dynamicPythonCode = dynamicPythonCode;
        this.startupTimeout = startupTimeout;
    }

    public boolean isDynamicPythonCode() {
        return dynamicPythonCode;
    }

    /**
     * Executable relative to AMUSE distribution root, absolute path in filesystem, or name of dynamic python code script
     */
    public String getExecutable() {
        return executable;
    }

    public String getWorkerDir() {
        return workerDir;
    }

    public int getNrOfWorkers() {
        return nrOfWorkers;
    }

    public int getNrOfThreads() {
        return nrOfThreads;
    }

    @Override
    public int getNrOfSlots() {
        return getNrOfWorkers();
    }

    public int getStartupTimeout() {
        return startupTimeout;
    }

    @Override
    public String getType() {
        return "worker";
    }

    /* (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        return "WorkerJobDescription [executable=" + executable + ", workerDir=" + workerDir + ", nrOfWorkers=" + nrOfWorkers
                + ", nrOfThreads=" + nrOfThreads + ", dynamicPythonCode=" + dynamicPythonCode + ", startupTimeout="
                + startupTimeout + ", id=" + id + ", stdoutFile=" + stdoutFile + ", stderrFile=" + stderrFile + ", label="
                + label + "]";
    }

}
