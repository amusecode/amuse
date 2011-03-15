package ibis.amuse;

import java.io.File;
import java.util.Properties;
import java.util.UUID;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.deploy.Deploy;
import ibis.deploy.Job;

public class Deployment {

    private static final Logger logger = LoggerFactory
            .getLogger(Deployment.class);
   
    private Deploy deploy;
    
    public Deployment() throws Exception {
       deploy = new Deploy(new File("deploy"));
       
    }
    
    public String getServerAddress() throws Exception {
        return deploy.getServerAddress();
    }

    public Job deploy(String codeName, String hostname, UUID workerID) {
        logger.info("Please deploy worker with ID " + workerID + " running "
                + codeName + " on host " + hostname);

        // TODO implement deployment
        return null;
    }

}
