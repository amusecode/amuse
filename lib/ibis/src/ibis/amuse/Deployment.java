package ibis.amuse;

import java.io.File;
import java.util.UUID;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.deploy.Application;
import ibis.deploy.Deploy;
import ibis.deploy.Experiment;
import ibis.deploy.Grid;
import ibis.deploy.Job;
import ibis.deploy.JobDescription;

public class Deployment {

    private static final Logger logger = LoggerFactory
            .getLogger(Deployment.class);

    private Deploy deploy;

    private final Grid grid;
    private final Experiment experiment;

    public Deployment() throws Exception {
        grid = new Grid(new File("deploy.grid"));
        experiment = new Experiment("amuse");

        deploy = new Deploy(new File("deploy"), true, true, 0, null, null, true);

    }

    public String getServerAddress() throws Exception {
        return deploy.getServerAddress();
    }

    public Job deploy(String codeName, String hostname, UUID workerID)
            throws Exception {
        logger.info("Deploying worker with ID " + workerID + " running \""
                + codeName + "\" on host " + hostname);


        JobDescription jobDescription = experiment.createNewJob(workerID
                .toString());

        if (hostname.equalsIgnoreCase("localhost")) {
            hostname = "local";
        }

        jobDescription.setClusterName(hostname);
        jobDescription.setProcessCount(1);
        jobDescription.setResourceCount(1);
        jobDescription.setRuntime(60);
        jobDescription.setApplicationName(codeName);

        Application application = jobDescription.getApplicationOverrides();

        application.setLibs(new File("deploy/lib-server"), new File("lib"));
        application.setMainClass("ibis.amuse.RemoteWorker");
        application.setMemorySize(1000);
        application.setLog4jFile(new File("log4j.properties"));

        application.setArguments("--code-name", codeName, "--worker-id",
                workerID.toString());

        Job result = deploy.submitJob(jobDescription, null, grid, null, null);

        result.waitUntilDeployed();

        return result;

    }
}
