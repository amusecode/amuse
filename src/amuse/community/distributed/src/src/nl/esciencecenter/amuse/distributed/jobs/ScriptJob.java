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

import ibis.ipl.Ibis;
import ibis.ipl.ReadMessage;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.util.FileTransfers;

/**
 * @author Niels Drost
 * 
 */
public class ScriptJob extends AmuseJob {

    private final ScriptJobDescription description;

    public ScriptJob(ScriptJobDescription description, Ibis ibis, JobSet jobManager) throws DistributedAmuseException {
        super(description, ibis, jobManager);
        this.description = description;
    }

    public ScriptJobDescription getDescription() {
        return description;
    }

    @Override
    void writeJobData(WriteMessage writeMessage) throws IOException {
        Path cwd = Paths.get("").toAbsolutePath();

        FileTransfers.writeDirectory(description.getScriptDir(), cwd, writeMessage);

        if (description.getInputDir() != null) {
            FileTransfers.writeDirectory(description.getInputDir(), cwd, writeMessage);
        }
    }

    @Override
    void readJobResult(ReadMessage readMessage) throws ClassNotFoundException, IOException {

        if (description.getOutputDir() != null) {
            Path cwd = Paths.get("").toAbsolutePath();
            FileTransfers.readDirectory(cwd, readMessage);
        }
    }

}
