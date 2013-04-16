package ibis.amuse;

import java.io.File;
import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.deploy.Application;
import ibis.deploy.ApplicationSet;
import ibis.deploy.Resource;
import ibis.deploy.Deploy;
import ibis.deploy.Experiment;
import ibis.deploy.Jungle;
import ibis.deploy.Job;
import ibis.deploy.JobDescription;
import ibis.deploy.Workspace;
import ibis.deploy.Deploy.HubPolicy;
import ibis.deploy.gui.GUI;
import ibis.deploy.gui.Mode;

public class Deployment {

    private static final Logger logger = LoggerFactory.getLogger(Deployment.class);

    public static final String[] logos = { "images/strw-logo-blue.png", "images/nova-logo.png" };

    private int nextMachine = 0;

    private final Deploy deploy;

    private final Jungle jungle;
    private final ApplicationSet applications;
    private final Experiment experiment;

    private final File amuseHome;
    private final File ibisDir;
    private final File logDir;

    private final int runtime;

    public Deployment(boolean verbose, boolean keepSandboxes, boolean useGui, boolean useHubs, File[] jungleFiles,
            String[] hubs, File logDir, int runtime) throws Exception {
        this.logDir = logDir;
        this.runtime = runtime;

        jungle = new Jungle();
        if (jungleFiles.length == 0) {
            logger.warn("No Jungle files specified, only local resource available");
        } else {
            for (File file : jungleFiles) {
                if (!file.exists()) {
                    throw new Exception("jungle description file " + file.getAbsolutePath() + " does not exist.");
                }
                if (!file.isFile() && !file.canRead()) {
                    throw new Exception("cannot read jungle description file " + file.getAbsolutePath());
                }

                jungle.load(file, false);

            }
        }

        if (logger.isInfoEnabled()) {
            Resource[] resources = jungle.getResources();
            String[] names = new String[resources.length];
            for (int i = 0; i < resources.length; i++) {
                names[i] = resources[i].getName();
            }

            logger.info("loaded " + jungleFiles.length + " files, resources available: " + Arrays.toString(names));
        }

        experiment = new Experiment("amuse");
        applications = new ApplicationSet();

        // location of AMUSE
        String amuseHomeProperty = System.getProperty("amuse.home");
        if (amuseHomeProperty == null) {
            throw new Exception("amuse.home property not specified");
        }

        amuseHome = new File(amuseHomeProperty).getAbsoluteFile();
        if (!amuseHome.exists()) {
            throw new Exception("amuse home (" + amuseHome + ") does not exist");
        }
        if (!amuseHome.isDirectory()) {
            throw new Exception("amuse home (" + amuseHome + ") is not a directory");
        }

        // location of Ibis Library in AMUSE
        ibisDir = new File(amuseHome, "lib/ibis");
        if (!ibisDir.exists()) {
            throw new Exception("ibis library dir (" + ibisDir + ") does not exist");
        }
        if (!ibisDir.isDirectory()) {
            throw new Exception("ibis library dir (" + ibisDir + ") is not a directory");
        }

        deploy = new Deploy(null, verbose, 0, null, null, true);
        deploy.setKeepSandboxes(keepSandboxes);
        deploy.setMonitoringEnabled(useGui);
        if (!useHubs) {
            deploy.setHubPolicy(HubPolicy.OFF);
        }

        for (String hub : hubs) {
            logger.info("Starting hub on " + hub);
            Resource hubResource = jungle.getResource(hub);

            if (hubResource == null) {
                throw new Exception("hub cannot be started on " + hub
                        + " as it is not specified in the resource description (jungle)");
            }

            deploy.getHub(hubResource, true, null);
        }

        // set location of wrapper scripts to ibis lib dir, IF it does not exist
        // in CWD
        for (Resource resource : jungle.getResources()) {
            if (resource.getJobWrapperScript() != null && !resource.getJobWrapperScript().exists()) {
                // wrapper script not in CWD, try some alternative locations

                File deployHomeLocation = new File(deploy.getHome(), resource.getJobWrapperScript().getName());
                File ibisDirLocation = new File(ibisDir, resource.getJobWrapperScript().getName());

                if (!deployHomeLocation.exists() && !ibisDirLocation.exists()) {
                    throw new Exception(resource.getJobWrapperScript().getName() + " not found in "
                            + resource.getJobWrapperScript().getAbsolutePath() + ", " + ibisDirLocation + ", or "
                            + deployHomeLocation + ", Specified in resource \"" + resource.getName() + "\"");
                } else if (ibisDirLocation.exists()) {
                    resource.setJobWrapperScript(ibisDirLocation);
                } else if (deployHomeLocation.exists()) {
                    resource.setJobWrapperScript(deployHomeLocation);
                }

            }
        }

        Workspace workspace = new Workspace(jungle, applications, experiment);

        if (useGui) {
            new GUI(deploy, workspace, Mode.MONITORING_ONLY, true, logos);
        }
    }

    /**
     * Return next RoundRobin machine. Skips local unless there is only one
     * machine
     * 
     * @return the next machine, selected round-robin
     */
    private synchronized Resource getNextMachine() {
        Resource[] resources = jungle.getResources();

        if (resources.length == 0) {
            return null;
        } else if (resources.length == 1) {
            return resources[0];
        }

        Resource result;

        do {
            if (nextMachine < 1 || nextMachine >= resources.length) {
                nextMachine = 0;
            }
            result = resources[nextMachine++];
        } while (result.getName().equals("local"));

        return result;
    }

    /**
     * Return a random machine. Will not return local local unless there is only
     * one machine
     * 
     * @return a resource selected at random
     */
    private synchronized Resource getRandomMachine() {
        Resource[] resources = jungle.getResources();

        if (resources.length == 0) {
            return null;
        } else if (resources.length == 1) {
            return resources[0];
        }

        Resource result;
        do {
            int selected = (int) (Math.random() * resources.length);
            result = resources[selected];
        } while (result.getName().equals("local"));

        return result;
    }

    public String getServerAddress() throws Exception {
        return deploy.getServerAddress();
    }

    public synchronized Job deploy(String codeName, String codeDir, String resourceName, String stdoutFile,
            String stderrFile, String workerID, int nrOfWorkers, int nrOfNodes, int nrOfThreads, boolean copyWorkerCode, String nodeLabel)
            throws Exception {
        Resource resource = null;
        logger.info("Deploying worker \"" + workerID + "\" on host " + resourceName + " with " + nrOfWorkers
                + " workers on " + nrOfNodes + " nodes with label " + nodeLabel);

        if (resourceName.equalsIgnoreCase("localhost")) {
            resourceName = "local";
        }
        if (resourceName.equals("random")) {
            resource = getRandomMachine();

            Resource[] resources = jungle.getResources();
            int selected = (int) (Math.random() * resources.length);
            resource = resources[selected];

            resourceName = resource.getName();
            logger.info("Randomly selected resource " + resourceName);
        } else if (resourceName.equals("roundrobin")) {
            resource = getNextMachine();

            resourceName = resource.getName();
            logger.info("Round-Robin selected resource " + resourceName);

        } else {
            resource = jungle.getResource(resourceName);
        }

        if (resource == null) {
            throw new Exception("Resource \"" + resourceName + "\" not found in jungle description");
        }

        String remoteAmuseHome = resource.getProperties().getProperty("amuse.home");

        if (resourceName.equals("local")) {
            remoteAmuseHome = new File(System.getProperty("amuse.home")).getAbsolutePath();
        }

        if (remoteAmuseHome == null) {
            throw new Exception("amuse.home property not set for resource \"" + resourceName
                    + "\" in jungle description");
        }

        // get or create Application for worker
        Application application = applications.getApplication(workerID);

        if (application == null) {
            application = new Application(workerID);
            applications.addApplication(application);

            File serverLib = new File(deploy.getHome(), "lib-server");

            if (!serverLib.isDirectory()) {
                throw new Exception(serverLib + " does not exist");
            }

            File ibisLibDir = new File(ibisDir, "lib");

            if (!ibisLibDir.isDirectory()) {
                throw new Exception(serverLib + " does not exist");
            }

            application.setLibs(serverLib, ibisLibDir);

            // application.addInputFile(new
            // File("libibis-amuse-bhtree_worker.so"));
            application.setMainClass("ibis.amuse.CodeProxy");
            application.setMemorySize(1000);
            application.setLog4jFile(new File(ibisDir, "log4j.properties"));

            application.setSystemProperty("ibis.managementclient", "true");
            application.setSystemProperty("ibis.bytescount", "true");

            if (copyWorkerCode) {
                File codeDirFile = new File(codeDir);
                application.addInputFile(new File(codeDirFile, codeName));

                if (!codeDirFile.isDirectory()) {
                    throw new Exception("cannot copy codedir " + codeDir
                            + " as it is not a directory, or does not exist");
                }

                for (File file : codeDirFile.listFiles()) {
                    application.addInputFile(file);
                }
            }
        }

        // create job description
        JobDescription jobDescription = new JobDescription(workerID);
        experiment.addJob(jobDescription);

        jobDescription.getResource().setName(resourceName);

        jobDescription.setProcessCount(nrOfNodes);
        jobDescription.setResourceCount(nrOfNodes);
        jobDescription.setRuntime(runtime);
        jobDescription.getApplication().setName(codeName);
        jobDescription.setPoolName("amuse");
        jobDescription.setStdoutFile(new File(logDir, workerID + ".out.txt"));
        jobDescription.setStderrFile(new File(logDir, workerID + ".err.txt"));

        // jobDescription.setPoolSize(nrOfNodes);

        jobDescription.getApplication().addOutputFile(new File("output"));

        String absCodeDir = remoteAmuseHome + "/" + codeDir;

        jobDescription.getApplication().setSystemProperty("java.library.path", absCodeDir);

        jobDescription.getApplication().setArguments("--code-name", codeName, "--worker-id", workerID, "--amuse-home",
                remoteAmuseHome, "--code-dir", codeDir, "--number-of-processes", Integer.toString(nrOfWorkers),
                "--number-of-nodes", Integer.toString(nrOfNodes), "--number-of-threads", Integer.toString(nrOfThreads),"--stdout", stdoutFile, "--stderr", stderrFile,
                "--copy-worker-code", "" + copyWorkerCode);

        Job result = deploy.submitJob(jobDescription, application, resource, null, null);

        return result;

    }

}
