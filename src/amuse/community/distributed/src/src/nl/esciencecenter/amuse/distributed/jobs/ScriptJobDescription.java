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
public class ScriptJobDescription extends AmuseJobDescription implements Serializable {

    private static final long serialVersionUID = 1L;

    private final String scriptName;
    private final String arguments;
    private final String scriptDir;
    private final String inputDir;
    private final String outputDir;

    public ScriptJobDescription(String stdoutFile, String stderrFile, String nodeLabel, String scriptName, String arguments,
            String scriptDir, String inputDir, String outputDir) {
        super(stdoutFile, stderrFile, nodeLabel);

        this.scriptName = scriptName;
        this.arguments = arguments;
        this.scriptDir = scriptDir;
        this.inputDir = inputDir;
        this.outputDir = outputDir;
    }

    @Override
    public int getNrOfSlots() {
        return 1;
    }

    @Override
    public String getType() {
        return "script";
    }

    public String getScriptName() {
        return scriptName;
    }

    public String getArguments() {
        return arguments;
    }

    public String getScriptDir() {
        return scriptDir;
    }

    public String getInputDir() {
        return inputDir;
    }

    public String getOutputDir() {
        return outputDir;
    }

    @Override
    public String toString() {
        return "ScriptJobDescription [scriptName=" + scriptName + ", arguments=" + arguments + ", scriptDir=" + scriptDir
                + ", inputDir=" + inputDir + ", outputDir=" + outputDir + ", ID=" + id + ", stdoutFile=" + stdoutFile
                + ", stderrFile=" + stderrFile + ", label=" + label + "]";
    }

}
