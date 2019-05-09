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
package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.Server;

import java.io.InputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.adaptors.schedulers.ssh.SshSchedulerAdaptor;
import nl.esciencecenter.xenon.credentials.Credential;
import nl.esciencecenter.xenon.credentials.DefaultCredential;
import nl.esciencecenter.xenon.filesystems.FileSystem;
import nl.esciencecenter.xenon.utils.LocalFileSystemUtils;
import nl.esciencecenter.xenon.filesystems.Path;
import nl.esciencecenter.xenon.schedulers.Scheduler;
import nl.esciencecenter.xenon.adaptors.schedulers.slurm.SlurmSchedulerAdaptor;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Manages a resource. Copies files, possibly starts a hub, etc.
 * 
 * @author Niels Drost
 * 
 */
public class ResourceManager {

    public static final String WHITESPACE_REGEX = ";";

    public static final String EQUALS_REGEX = "\\s*=\\s*";

    private static final Logger logger = LoggerFactory.getLogger(ResourceManager.class);

    private static int nextID = 0;

    private static int getNextID() {
        return nextID++;
    }

    private final int id;
    private final String name;
    private final String username;
    private final String location;
    private final String gateway;
    private final String amuseDir;
    private final String tmpDir;
    private final String schedulerType;
    private final String hubQueueName;
    private final int hubTimeMinutes;

    private final Map<String, String> properties = new HashMap<String, String>();
    private final Credential credential;
    
    private final AmuseConfiguration configuration;

    private final boolean startHub;

    private final Scheduler scheduler;
    private final Path home;
    private final FileSystem filesystem;

    private final String keyfile;

    private final Hub hub;

    private static void waitUntilHubStarted(Server iplServer, String hubAddress, String name) throws DistributedAmuseException {
        logger.info("waiting for new remote hub on {} at {} to connect to the local hub", name, hubAddress);
        for (int i = 0; i < 40; i++) {
            String[] knownHubAddresses = iplServer.getHubs();
            logger.debug("ipl hub addresses now " + Arrays.toString(iplServer.getHubs()));
            for (String knownHub : knownHubAddresses) {
                if (knownHub.equals(hubAddress)) {
                    logger.info("new hub at {} now connected to local hub", name);
                    return;
                }
            }
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                //IGNORE
            }
        }
        throw new DistributedAmuseException("Local and new remote Hub at " + name + " not able to communicate");
    }

    public ResourceManager(String name, String username, String location, String gateway, String amuseDir, String tmpDir, String schedulerType, String hubQueueName, int hubTimeMinutes, boolean startHub,
            Server iplServer) throws DistributedAmuseException {
        this.id = getNextID();
        this.name = name;
        this.location = location;
        this.gateway = gateway;
        this.amuseDir = amuseDir;
        this.tmpDir = tmpDir;
        this.schedulerType = schedulerType;
        this.hubQueueName = hubQueueName;
        this.hubTimeMinutes = hubTimeMinutes;
        this.username = username;
        this.keyfile = "id_rsa"; // hardcoded and non functional atm
        
        //local resources _never_ have a hub
        this.startHub = (schedulerType.equals("local")) ? false : startHub;
        
        if (gateway != null && !gateway.isEmpty()) {
            properties.put(SshSchedulerAdaptor.GATEWAY, gateway);
        }
        //~ if(getSchedulerType().equals("slurm")) {
            //~ properties.put(SlurmSchedulerAdaptor.IGNORE_VERSION_PROPERTY, "true");
        //~ }

        if(username != null) {
            credential = new DefaultCredential(username);
        } else {
            credential = new DefaultCredential();
        } 

        filesystem = _getFileSystem();
        
        home = _getHome();
        logger.info("found home of resource {} to be {}", name, home);        

        this.configuration = downloadConfiguration(filesystem);

        if (!configuration.isJavaEnabled()) {
            throw new DistributedAmuseException("Resource " + name
                    + " not suitable as target for distributed AMUSE, java not enabled in configuration");
        }

        scheduler = createScheduler();

        if (this.startHub) {
            this.hub = new Hub(this, this.configuration, iplServer.getHubs());
            iplServer.addHubs(this.hub.getAddress());

            String hubAddress = this.hub.getAddress();
            logger.debug("just added new hub " + hubAddress);

            waitUntilHubStarted(iplServer, hubAddress, name);
        } else {
            this.hub = null;
        }
        logger.info("Created new resource {}", this);

    }

    private Scheduler createScheduler() throws DistributedAmuseException {
        try {
            if (isLocal()) {
                return Scheduler.create("local");
            }

            return Scheduler.create(getSchedulerType(), getLocation(), getCredential(), getProperties());
        } catch (XenonException e) {
            throw new DistributedAmuseException("cannot create scheduler connection for resource " + this.name, e);
        }
    }

    private Path _getHome() throws DistributedAmuseException {
        return getFileSystem().getWorkingDirectory();
    }

    private FileSystem _getFileSystem() throws DistributedAmuseException {
        try {
            if (isLocal()) {
                FileSystem filesystem = LocalFileSystemUtils.getLocalFileSystems()[0];
                filesystem.setWorkingDirectory(new Path(System.getProperty("user.home")));
                return filesystem;
            }


            logger.info("trying with {} {} {}", username, getLocation(), getProperties());

            return FileSystem.create("sftp", getLocation(), getCredential(), getProperties());
        } catch (XenonException e) {
            throw new DistributedAmuseException("cannot open filesystem for resource " + this.name, e);
        }
    }

    private AmuseConfiguration downloadConfiguration(FileSystem filesystem) throws DistributedAmuseException {
        try {
            Path amuseHome;
            if (this.amuseDir.startsWith("/")) {
                amuseHome = new Path(this.amuseDir);
            } else {
                Path userHome = filesystem.getWorkingDirectory();

                amuseHome = userHome.resolve(this.amuseDir);
            }
            
            Path amuseConfig = amuseHome.resolve("config.mk");
            if(!filesystem.exists(amuseConfig))
                amuseConfig = amuseHome.resolve("../../../../share/amuse/config.mk");
            if(!filesystem.exists(amuseConfig)) 
                amuseConfig = amuseHome.resolve("../../../../../share/amuse/config.mk");
            if(!filesystem.exists(amuseConfig))
                amuseConfig = amuseHome.resolve("../../config.mk");
            if(!filesystem.exists(amuseConfig))
                throw new DistributedAmuseException("cannot find config file config.mk from " + amuseHome);

            logger.debug("Downloading amuse config for " + getName() + " from " + amuseConfig);

            InputStream in = filesystem.readFromFile(amuseConfig);

            return new AmuseConfiguration(amuseHome.toAbsolutePath().toString(), in);
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot download configuration file for resource " + this.name, e);
        }
    }

    public int getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public String getLocation() {
        return location;
    }

    public String getGateway() {
        return gateway;
    }

    public String getAmuseDir() {
        return amuseDir;
    }
    
    public String getTmpDir() {
        return tmpDir;
    }


    public String getSchedulerType() {
        return schedulerType;
    }

    public AmuseConfiguration getConfiguration() {
        return configuration;
    }

    @Override
    public int hashCode() {
        return new Integer(id).hashCode();
    }

    @Override
    public boolean equals(Object other) {
        if (other == null) {
            return false;
        }

        if (!(other instanceof ResourceManager)) {
            return false;
        }

        return id == ((ResourceManager) other).id;
    }

    public void stop() {
        logger.debug("Stopping resource {}", this);
        if (hub != null) {
            hub.stopHub();
        }
        try {
            scheduler.close();
        } catch (XenonException e) {
            logger.warn("Error while closing scheduler for " + this, e);
        }
        try {
            filesystem.close();
        } catch (XenonException e) {
            logger.warn("Error while closing filesystem for " + this, e);
        }

    }

    public Hub getHub() {
        return hub;
    }

    public boolean hasHub() {
        return hub != null;
    }

    public boolean isLocal() {
        return location == null || location.equals("localhost") || location.equals("local");
    }

    @Override
    public String toString() {
        return "Resource [id=" + id + ", name=" + name + ", username=" + username + ", location=" + location + ", amuseDir=" + amuseDir + ", schedulerType="
                + schedulerType + ", configuration=" + configuration + ", startHub=" + startHub + ", hub=" + hub + "]";
    }

    public Map<String, String> getStatusMap() throws DistributedAmuseException {
        Map<String, String> result = new LinkedHashMap<String, String>();

        result.put("ID", Integer.toString(id));
        result.put("Name", name);
        result.put("Username", username);
        result.put("Location", location);
        result.put("Gateway", gateway);
        result.put("Amuse dir", amuseDir);
        result.put("Scheduler type", schedulerType);

        result.put("Java path", configuration.getJava());
        result.put("MPI enabled", Boolean.toString(configuration.isMPIEnabled()));
        result.put("Mpiexec", configuration.getMpiexec());

        return result;
    }

    public Path getHome() {
        return home;
    }

    public Scheduler getScheduler() {
        return scheduler;
    }

    public int getHubTimeMinutes() {
        return this.hubTimeMinutes;
    }

    public String getHubQueueName() {
        return this.hubQueueName;
    }

    public FileSystem getFileSystem() {
        return filesystem;
    }

    public Map<String, String> getProperties() {
        return properties;
    }

    public Credential getCredential() {
        return credential;
    }

    public String getKeyfile() {
        return keyfile;
    }

}
