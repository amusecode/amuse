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
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.lang.reflect.Field;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Arrays;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.jobs.AmuseJobDescription;
import nl.esciencecenter.amuse.distributed.workers.OutputForwarder;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public abstract class JobRunner extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(JobRunner.class);

    protected final AmuseJobDescription description;
    protected final AmuseConfiguration amuseConfiguration;
    protected final Ibis ibis;
    protected final ReceivePortIdentifier resultPort;
    protected final Path tmpDir;
    protected final Path sandbox;
    protected final boolean debug;

    private Process process = null;
    private OutputForwarder out;
    private OutputForwarder err;

    private Exception error = null;

    public JobRunner(AmuseJobDescription description, AmuseConfiguration amuseConfiguration, ReceivePortIdentifier resultPort,
            Ibis ibis, Path tmpDir,boolean debug) throws IOException {
        this.description = description;
        this.amuseConfiguration = amuseConfiguration;
        this.resultPort = resultPort;
        this.ibis = ibis;
        this.tmpDir = tmpDir;
        this.debug = debug;

        this.sandbox = tmpDir.resolve("job-" + description.getID());
        Files.createDirectory(this.sandbox);
    }

    protected synchronized void setError(Exception error) {
        if (this.error == null) {
            this.error = error;
        } else {
            logger.warn("Masked error in running job: " + error);
        }
    }

    protected synchronized Exception getError() {
        return error;
    }

    private synchronized Process getProcess() {
        return process;
    }

    synchronized void startProcess(ProcessBuilder builder) throws IOException {
        logger.debug("Running process in cwd: " + builder.directory());

        process = builder.start();

        //attach streams
        out = new OutputForwarder(process.getInputStream(), description.getStdoutFile(), ibis);
        err = new OutputForwarder(process.getErrorStream(), description.getStderrFile(), ibis);

    }

    int getExitCode() {
        return getProcess().exitValue();
    }

    void waitForProcess() {
        Process process = getProcess();

        try {
            int result = process.waitFor();

            if (result != 0) {
                setError(new Exception("Process ended with non-zero exit code: " + result));
            }
        } catch (InterruptedException e) {
            //IGNORE
        }

        out.waitFor(1000);
        err.waitFor(1000);
    }

    private void nativeKill() {
        Process process = getProcess();

        if (process == null) {
            return;
        }

        try {
            Field f = process.getClass().getDeclaredField("pid");
            f.setAccessible(true);

            Object pid = f.get(process);

            ProcessBuilder builder = new ProcessBuilder("/bin/sh", "-c", "kill -9 " + pid.toString());

            builder.redirectError(Redirect.INHERIT);
            //builder.redirectInput();
            builder.redirectOutput(Redirect.INHERIT);

            logger.info("Killing process using command: " + Arrays.toString(builder.command().toArray()));

            Process killProcess = builder.start();

            killProcess.getOutputStream().close();

            int exitcode = killProcess.waitFor();

            logger.info("native kill done, kill exit code is " + exitcode);

        } catch (Throwable t) {
            logger.error("Error on (forcibly) killing process", t);
        }
    }

    protected void deleteSandbox() {
        if (!Files.exists(sandbox)) {
            logger.debug("Sandbox does not exist, not cleaning up");
            return;
        }
        
        if (debug) {
            logger.debug("Not cleaning up sandbox because we are in debug mode");
            return;
        }

        logger.debug("Cleaning up Sandbox");

        try {
            Files.walkFileTree(sandbox, new SimpleFileVisitor<Path>() {
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                    Files.delete(file);
                    return FileVisitResult.CONTINUE;
                }

                @Override
                public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                    Files.delete(dir);
                    return FileVisitResult.CONTINUE;
                }

            });
        } catch (IOException e) {
            logger.error("Error while cleaning up sandbox at: " + sandbox, e);
        }
    }

    //small utility to figure out if the process is still running.
    protected boolean hasEnded() {
        Process process = getProcess();

        try {
            int exitcode = process.exitValue();
            //we only end up here if the process is done

            logger.info("Process ended with value " + exitcode);
            return true;
        } catch (IllegalThreadStateException e) {
            //we got this exception as the process is not done yet
            return false;
        }
    }

    protected void sendResult() {
        logger.debug("Job done. Sending result to main amuse node.");

        //send result message to job
        try {
            SendPort sendPort = ibis.createSendPort(DistributedAmuse.MANY_TO_ONE_PORT_TYPE);

            sendPort.connect(resultPort);

            WriteMessage message = sendPort.newMessage();

            message.writeObject(getError());

            writeResultData(message);

            message.finish();

            sendPort.close();

            logger.debug("result sent.");
        } catch (IOException e) {
            logger.error("Failed to report status to main node", e);
        }

    }

    /**
     * This function will cancel the job by killing the process running the job. Since it is a common case to cancel a job when a
     * normally finishing job is cleaned up, we give the process a short while to finish on its own.
     */
    public synchronized void cancel() {
        try {

            if (process == null) {
                return;
            }

            //wait a little while for the process (or actually the streams connected to the process) to finish by itself
            out.waitFor(50);
            err.waitFor(50);

            //check if the process is finished
            if (hasEnded()) {
                return;
            }
            ;

            process.destroy();

            //wait a little while for the process (or actually the streams connected to the process) to finish.
            out.waitFor(50);
            err.waitFor(50);

            if (hasEnded()) {
                return;
            }
            ;

            logger.error("Process not ended after process.destroy()! Performing native kill");
            nativeKill();
        } finally {
            deleteSandbox();
        }
    }

    abstract void writeResultData(WriteMessage message) throws IOException;

}
