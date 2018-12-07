package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.Server;
import ibis.ipl.server.ServerProperties;

import java.io.File;
import java.util.ArrayList;
import java.util.Properties;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Collection of resources potentially available for starting reservations on.
 * 
 * @author Niels Drost
 * 
 */
public class ResourceSet {

    private static final Logger logger = LoggerFactory.getLogger(ResourceSet.class);

    private final Server iplServer;

    private final ArrayList<ResourceManager> resources;

    private final boolean startHubs;
    
    private final String amuseRootDir;

    public ResourceSet(String amuseRootDir, boolean startHubs) throws DistributedAmuseException {
        resources = new ArrayList<ResourceManager>();
        this.startHubs = startHubs;
        this.amuseRootDir = amuseRootDir;

        try {
            Properties properties = new Properties();
            //use a random free port.
            properties.put(ServerProperties.PORT, "0");
            this.iplServer = new Server(properties);
        } catch (Exception e) {
            throw new DistributedAmuseException("could not create IPL server", e);
        }

        //add local resource by default

        logger.debug("local amuse dir = " + amuseRootDir);
        //newResource("local", null, null, amuseRootDir, "local");
    }

    public synchronized ResourceManager newResource(String name, String location, String gateway, String amuseDir, String tmpDir,
            String schedulerType, String queueName, int timeMinutes) throws DistributedAmuseException {
        logger.debug("creating new resource: name = " + name + " location = " + location + " scheduler type = " + schedulerType
                + " amuse dir = " + amuseDir);

        for (ResourceManager resource : resources) {
            if (resource.getName().equals(name)) {
                throw new DistributedAmuseException("Resource " + name + " already exists");
            }
        }

        String gatewayLocation = null;
        if (gateway != null && !gateway.isEmpty()) {
            gatewayLocation = getResource(gateway).getLocation();
        }

        String username = null;
        String hostname = location;
        if(location != null && location.indexOf("@") != -1 ) {
            username = location.substring(0, location.indexOf("@"));
            hostname = location.substring(location.indexOf("@")+1); 
          }

        ResourceManager result = new ResourceManager(name, username, hostname, gatewayLocation, amuseDir, tmpDir, 
                                      schedulerType, queueName, timeMinutes, this.startHubs, iplServer);

        resources.add(result);

        return result;
    }

    public synchronized ResourceManager getResource(int resourceID) throws DistributedAmuseException {
        for (ResourceManager resource : resources) {
            if (resource.getId() == resourceID) {
                return resource;
            }
        }
        throw new DistributedAmuseException("Resource with ID " + resourceID + " not found");
    }

    public synchronized ResourceManager getResource(String name) throws DistributedAmuseException {
        for (ResourceManager resource : resources) {
            if (resource.getName().equals(name)) {
                return resource;
            }
        }
        throw new DistributedAmuseException("Resource with name " + name + " not found");
    }

    public synchronized int getResourceCount() {
        return resources.size();
    }

    public synchronized ResourceManager[] getResources() {
        return resources.toArray(new ResourceManager[resources.size()]);
    }

    public synchronized void deleteResource(ResourceManager resource) throws DistributedAmuseException {
        for (int i = 0; i < resources.size(); i++) {
            if (resource.getId() == resources.get(i).getId()) {
                resource.stop();
                resources.remove(i);
                return;
            }
        }
        throw new DistributedAmuseException("Resource " + resource.getId() + " not found");
    }

    public String getIplServerAddress() {
        return iplServer.getAddress();
    }

    public String[] getHubAddresses() {
        return iplServer.getHubs();
    }

    public void endRegistry() {
        iplServer.getRegistryService().end(60000);
    }

    public synchronized void end() {
        logger.debug("waiting for all ipl services to end");
        iplServer.end(60000);
        logger.debug("services ended");

        for (ResourceManager resource : resources) {
            logger.debug("ending resource {}", resource);
            resource.stop();
        }
    }

}
