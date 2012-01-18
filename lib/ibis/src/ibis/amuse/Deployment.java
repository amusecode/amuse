package ibis.amuse;

import java.io.File;

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

    private final Deploy deploy;

    private final Jungle jungle;
    private final ApplicationSet applications;
    private final Experiment experiment;

    private final File amuseHome;
    private final File ibisDir;

    public Deployment(boolean verbose, boolean keepSandboxes, boolean useGui, boolean useHubs, String... logos)
            throws Exception {
        jungle = new Jungle(new File("deploy.jungle"), true);

        experiment = new Experiment("amuse");
        applications = new ApplicationSet();

        // location of AMUSE
        String amuseHomeProperty = System.getProperty("amuse.home");
        if (amuseHomeProperty == null) {
            throw new Exception("amuse.home propertyy not specified");
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

        deploy = new Deploy(null, verbose, keepSandboxes, useGui, true, 0, null, null, true);
        if (!useHubs) {
            deploy.setHubPolicy(HubPolicy.OFF);
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
                            + resource.getJobWrapperScript().getAbsolutePath() + ", " + ibisDirLocation + ", or " + deployHomeLocation
                            + ", Specified in resource \"" + resource.getName() + "\" of deploy.jungle");
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

    public String getServerAddress() throws Exception {
        return deploy.getServerAddress();
    }

    public Job deploy(String codeName, String codeDir, String resourceName, String workerID, int nrOfWorkers,
            int nrOfNodes) throws Exception {
        logger.info("Deploying worker \"" + workerID + "\" running \"" + codeName + "\" on host " + resourceName
                + " with " + nrOfWorkers + " workers on " + nrOfNodes + " nodes");

        if (resourceName.equalsIgnoreCase("localhost")) {
            resourceName = "local";
        }
        Resource resource = jungle.getResource(resourceName);

        if (resource == null) {
            throw new Exception("Resource \"" + resourceName
                    + "\"not found in jungle description file \"deploy.jungle\"");
        }

        String remoteAmuseHome = resource.getProperties().getProperty("amuse.home");

        if (resourceName.equals("local")) {
            remoteAmuseHome = new File(System.getProperty("amuse.home")).getAbsolutePath();
        }

        if (remoteAmuseHome == null) {
            throw new Exception("amuse.home property not set for resource \"" + resourceName
                    + "\" in jungle description file deploy.jungle");
        }

        String mpirun = resource.getProperties().getProperty("mpirun");

        if (mpirun == null) {
            if (!resourceName.equals("local")) {
                logger.warn("mpirun property not set for resource \"" + resourceName
                        + "\" in jungle description file deploy.jungle, using default (mpirun)");
            }
            mpirun = "mpirun";
        }

        // get or create Application for worker
        Application application = applications.getApplication(codeName);

        if (application == null) {
            application = new Application(codeName);
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
            application.setMainClass("ibis.amuse.CodeInterface");
            application.setMemorySize(1000);
            // application.setLog4jFile(new File("log4j.properties"));

            application.setSystemProperty("ibis.managementclient", "true");
            application.setSystemProperty("ibis.bytescount", "true");
        }

        // create job description
        JobDescription jobDescription = new JobDescription(workerID);
        experiment.addJob(jobDescription);

        jobDescription.getResource().setName(resourceName);

        jobDescription.setProcessCount(nrOfNodes);
        jobDescription.setResourceCount(nrOfNodes);
        jobDescription.setRuntime(60);
        jobDescription.getApplication().setName(codeName);
        jobDescription.setPoolName("amuse");

        jobDescription.getApplication().addOutputFile(new File("output"));

        String absCodeDir = remoteAmuseHome + "/" + codeDir;

        jobDescription.getApplication().setSystemProperty("java.library.path", absCodeDir);

        jobDescription.getApplication().setArguments("--code-name", codeName, "--worker-id", workerID, "--amuse-home",
                remoteAmuseHome, "--code-dir", codeDir, "--number-of-workers", Integer.toString(nrOfWorkers),
                "--number-of-nodes", Integer.toString(nrOfNodes), "--mpirun", mpirun);

        Job result = deploy.submitJob(jobDescription, application, resource, null, null);

        result.waitUntilDeployed();

        return result;

    }
}
