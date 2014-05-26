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
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJobDescription;
import nl.esciencecenter.amuse.distributed.jobs.ScriptJobDescription;
import nl.esciencecenter.amuse.distributed.util.FileTransfers;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class ScriptJobRunner extends JobRunner {

    private static final Logger logger = LoggerFactory.getLogger(ScriptJobRunner.class);

    private final ScriptJobDescription description;
    
    public ScriptJobRunner(AmuseJobDescription description, AmuseConfiguration configuration, ReceivePortIdentifier resultPort,
            Ibis ibis, Path tmpDir, boolean debug, ReadMessage message) throws Exception {
        super(description, configuration, resultPort, ibis, tmpDir, debug);

        this.description = (ScriptJobDescription) description;

        FileTransfers.readDirectory(sandbox, message);

        if (this.description.getInputDir() != null) {
            FileTransfers.readDirectory(sandbox, message);
        }

        message.finish();

        Path outputDir = sandbox.resolve(this.description.getOutputDir());

        Files.createDirectories(outputDir);

        startProcess(buildProcessBuilder());

        //start a thread to send back result when the job is done
        setName("Script Job Runner for Job " + description.getID());
        setDaemon(true);
        start();
    }

    private ProcessBuilder buildProcessBuilder() throws IOException, DistributedAmuseException {
        ProcessBuilder builder = new ProcessBuilder();

        builder.directory(sandbox.toFile());

        if (amuseConfiguration.isMPIEnabled()) {
            // use mpiexec
            builder.command().add(amuseConfiguration.getMpiexec());
            builder.command().add("-n");
            builder.command().add("1");
        }
        
        Path scriptPath = sandbox.resolve(description.getScriptDir()).resolve(description.getScriptName());

        Path amuseScriptPath = amuseConfiguration.getAmuseHome().getAbsoluteFile().toPath().resolve("amuse.sh");

        
        builder.command().add(amuseScriptPath.toString());

        builder.command().add(scriptPath.toString());

        if (!description.getArguments().isEmpty()) {
            for (String argument : description.getArguments().split("\\s")) {
                builder.command().add(argument);
            }
        }

        //remove some environment variables, if present
        builder.environment().remove("DISPLAY");

        logger.info("starting script process, command = " + builder.command());

        return builder;
    }

    @Override
    public void run() {
        waitForProcess();
    
        sendResult();

        deleteSandbox();

    }

    @Override
    void writeResultData(WriteMessage writeMessage) throws IOException {
        if (description.getOutputDir() != null) {
            FileTransfers.writeDirectory(description.getOutputDir(), sandbox, writeMessage);
        }
    }

}
