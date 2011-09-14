package ibis.amuse;

import java.io.File;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.deploy.Application;
import ibis.deploy.ApplicationSet;
import ibis.deploy.Cluster;
import ibis.deploy.Deploy;
import ibis.deploy.Experiment;
import ibis.deploy.Grid;
import ibis.deploy.Job;
import ibis.deploy.JobDescription;
import ibis.deploy.Workspace;
import ibis.deploy.gui.GUI;
import ibis.deploy.gui.Mode;

public class Deployment {

    private static final Logger logger = LoggerFactory
            .getLogger(Deployment.class);

    private final Deploy deploy;

    private final GUI gui;

    private final Grid grid;
    private final ApplicationSet applications;
    private final Experiment experiment;

    public Deployment(boolean verbose, boolean useGui, String... logos)
            throws Exception {
        grid = new Grid(new File("deploy.grid"));
        experiment = new Experiment("amuse");
        applications = new ApplicationSet();

        Workspace workspace = new Workspace(grid, applications, experiment);

        deploy = new Deploy(new File("deploy"), verbose, false, 0, null, null,
                true);

        if (useGui) {
            gui = new GUI(deploy, workspace, Mode.MONITOR, logos);
        } else {
            gui = null;
        }
    }

    public String getServerAddress() throws Exception {
        return deploy.getServerAddress();
    }

    public Job deploy(String codeName, String codeDir, String clusterName,
            String workerID, int nrOfWorkers, int nrOfNodes) throws Exception {
        logger.info("Deploying worker \"" + workerID + "\" running \""
                + codeName + "\" on host " + clusterName + " with " + nrOfWorkers + " workers on " + nrOfNodes + " nodes");

        if (clusterName.equalsIgnoreCase("localhost")) {
            clusterName = "local";
        }
        Cluster cluster = grid.getCluster(clusterName);

        if (cluster == null) {
            throw new Exception("Cluster \"" + clusterName
                    + "\"not found in grid description file \"deploy.grid\"");
        }

        String amuseHome = cluster.getProperties().getProperty("amuse.home");

        if (clusterName.equals("local")) {
            amuseHome = new File("../..").getAbsolutePath().toString();
        }

        if (amuseHome == null) {
            throw new Exception("amuse.home property not set for cluster \""
                    + clusterName + "\" in grid description file deploy.grid");
        }

        String mpirun = cluster.getProperties().getProperty("mpirun");

        if (mpirun == null) {
            if (!clusterName.equals("local")) {
                logger.warn("mpirun property not set for cluster \""
                        + clusterName
                        + "\" in grid description file deploy.grid, using default (mpirun)");
            }
            mpirun = "mpirun";
        }

        // get or create Application for worker
        Application application = applications.getApplication(codeName);

        if (application == null) {
            application = new Application(codeName);
            applications.addApplication(application);

            application.setLibs(new File("deploy/lib-server"), new File("lib"));

            // application.addInputFile(new
            // File("libibis-amuse-bhtree_worker.so"));
            application.setMainClass("ibis.amuse.CodeInterface");
            application.setMemorySize(1000);
            application.setLog4jFile(new File("log4j.properties"));

            application.setSystemProperty("ibis.managementclient", "true");
            application.setSystemProperty("ibis.bytescount", "true");
        }

        // create job description
        JobDescription jobDescription = new JobDescription(workerID);
        experiment.addJob(jobDescription);

        jobDescription.getCluster().setName(clusterName);

        jobDescription.setProcessCount(nrOfNodes);
        jobDescription.setResourceCount(nrOfNodes);
        jobDescription.setRuntime(60);
        jobDescription.getApplication().setName(codeName);
        jobDescription.setPoolName("amuse");

        jobDescription.getApplication().addOutputFile(new File("output"));
        
        String absCodeDir = amuseHome + "/" + codeDir;

        jobDescription.getApplication().setSystemProperty("java.library.path",
                absCodeDir);

        jobDescription.getApplication().setArguments("--code-name", codeName,
                "--worker-id", workerID, "--amuse-home", amuseHome,
                "--code-dir", codeDir, "--number-of-workers", 
                Integer.toString(nrOfWorkers),"--number-of-nodes",
                Integer.toString(nrOfNodes), "--mpirun", mpirun);

        Job result = deploy.submitJob(jobDescription, application, cluster,
                null, null);

        result.waitUntilDeployed();

        return result;

    }
}
