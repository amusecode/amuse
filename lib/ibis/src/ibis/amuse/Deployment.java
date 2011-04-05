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

    public Deployment(boolean verbose, boolean useGui) throws Exception {
        grid = new Grid(new File("deploy.grid"));
        experiment = new Experiment("amuse");
        applications = new ApplicationSet();

        Workspace workspace = new Workspace(grid, applications, experiment);

        deploy = new Deploy(new File("deploy"), verbose, false, 0, null, null,
                true);

        if (useGui) {
            gui = new GUI(deploy, workspace, Mode.MONITOR);
        } else {
            gui = null;
        }
    }

    public String getServerAddress() throws Exception {
        return deploy.getServerAddress();
    }

    public Job deploy(String codeName, String hostname, String workerID)
            throws Exception {
        logger.info("Deploying worker \"" + workerID + "\" running \""
                + codeName + "\" on host " + hostname);

        if (hostname.equalsIgnoreCase("localhost")) {
            hostname = "local";
        }
        Cluster cluster = grid.getCluster(hostname);

        if (cluster == null) {
            throw new Exception("Cluster \"" + hostname
                    + "\"not found in grid description file \"deploy.grid\"");
        }

        // get or create Application for worker
        Application application = applications.getApplication(codeName);

        if (application == null) {
        	application = new Application(codeName);
        	applications.addApplication(application);

            application.setLibs(new File("deploy/lib-server"), new File("lib"));

            // application.addInputFile(new
            // File("libibis-amuse-bhtree_worker.so"));
            application.setMainClass("ibis.amuse.RemoteWorker");
            application.setMemorySize(1000);
            application.setLog4jFile(new File("log4j.properties"));
        }

        // create job description
        JobDescription jobDescription = new JobDescription(workerID);
        experiment.addJob(jobDescription);
        

        jobDescription.getCluster().setName(hostname);
        jobDescription.setProcessCount(1);
        jobDescription.setResourceCount(1);
        jobDescription.setRuntime(60);
        jobDescription.getApplication().setName(codeName);
        jobDescription.setPoolName("amuse");

        String amuseHome = cluster.getProperties().getProperty("amuse.home");

        if (amuseHome == null) {
            throw new Exception("amuse.home property not set for cluster \""
                    + hostname + "\" in grid description file deploy.grid");
        }

        jobDescription.getApplication().setSystemProperty(
                "java.library.path", amuseHome + "/src/amuse/community/bhtree" + ":" + amuseHome + "/lib/ibis/src/native");

//        jobDescription.getApplicationOverrides().setArguments("--code-name",
//                codeName, "--worker-id", workerID);
        
      jobDescription.getApplication().setArguments("--code-name",
      "/home/niels/workspace/amuse/src/amuse/community/bhtree/bhtree_worker", "--worker-id", workerID);
        

        Job result = deploy.submitJob(jobDescription, application, cluster, null, null);

        result.waitUntilDeployed();

        return result;

    }
}
