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
import nl.esciencecenter.xenon.Xenon;
import nl.esciencecenter.xenon.XenonException;
import nl.esciencecenter.xenon.adaptors.ssh.SshAdaptor;
import nl.esciencecenter.xenon.credentials.Credential;
import nl.esciencecenter.xenon.files.FileSystem;
import nl.esciencecenter.xenon.files.Path;
import nl.esciencecenter.xenon.files.RelativePath;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class Resource {

    private static final Logger logger = LoggerFactory.getLogger(Resource.class);

    private static int nextID = 0;

    private static int getNextID() {
        return nextID++;
    }

    private final int id;
    private final String name;
    private final String location;
    private final String gateway;
    private final String amuseDir;
    private final String schedulerType;

    private final AmuseConfiguration configuration;

    //null == auto
    private final Boolean startHub;

    private final Hub hub;

    private static void waitUntilHubStarted(Server iplServer, String hubAddress, String name) throws DistributedAmuseException {
        logger.info("waiting for new remote hub on {} to connect to the local hub", name);
        for (int i = 0; i < 40; i++) {
            String[] knownHubAddresses = iplServer.getHubs();
            logger.trace("ipl hub addresses now " + Arrays.toString(iplServer.getHubs()));
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

    public Resource(String name, String location, String gateway, String amuseDir, String schedulerType, Boolean startHub,
            Xenon xenon, Server iplServer) throws DistributedAmuseException {
        this.id = getNextID();
        this.name = name;
        this.location = location;
        this.gateway = gateway;
        this.amuseDir = amuseDir;
        this.schedulerType = schedulerType;
        this.startHub = startHub;

        this.configuration = downloadConfiguration(xenon);

        if (!configuration.isJavaEnabled()) {
            throw new DistributedAmuseException("Resource " + name
                    + " not suitable as target for distributed AMUSE, java not enabled in configuration");
        }

        if (mustStartHub()) {
            this.hub = new Hub(this, this.configuration, iplServer.getHubs(), xenon);
            iplServer.addHubs(this.hub.getAddress());

            String hubAddress = this.hub.getAddress();
            logger.debug("just added new hub " + hubAddress);

            waitUntilHubStarted(iplServer, hubAddress, name);
        } else {
            this.hub = null;
        }
        logger.info("Created new resource {}", this);
    }

    private FileSystem openFileSystem(Xenon xenon) throws DistributedAmuseException {
        try {

            FileSystem filesystem;

            if (this.name.equals("local")) {
                filesystem = xenon.files().newFileSystem("local", "/", null, null);
            } else {
                Credential credential = xenon.credentials().getDefaultCredential("ssh");

                Map<String, String> properties = new HashMap<String, String>();
                String gateway = getGateway();
                if (gateway != null && !gateway.isEmpty()) {
                    properties.put(SshAdaptor.GATEWAY, gateway);
                }

                filesystem = xenon.files().newFileSystem("ssh", location, credential, properties);
            }

            return filesystem;
        } catch (XenonException e) {
            throw new DistributedAmuseException("cannot open filesystem for resource " + this.name, e);
        }
    }

    private AmuseConfiguration downloadConfiguration(Xenon xenon) throws DistributedAmuseException {
        FileSystem filesystem = openFileSystem(xenon);

        try {
            RelativePath amuseConfig = new RelativePath(this.amuseDir + "/config.mk");

            Path path = xenon.files().newPath(filesystem, amuseConfig);

            try (InputStream in = xenon.files().newInputStream(path)) {
                return new AmuseConfiguration(this.amuseDir, in);
            }
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot download configuration file for resource " + this.name, e);
        } finally {
            try {
                xenon.files().close(filesystem);
            } catch (XenonException e) {
                //IGNORE
            }
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

        if (!(other instanceof Resource)) {
            return false;
        }

        return id == ((Resource) other).id;
    }

    public boolean mustStartHub() {
        if (startHub == null) {
            return getSchedulerType() != null && (getSchedulerType().toLowerCase() != "local");
        }
        return startHub;
    }

    public void stop() {
        logger.debug("Stopping resource {}", this);
        if (hub != null) {
            hub.stop();
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
        return "Resource [id=" + id + ", name=" + name + ", location=" + location + ", amuseDir=" + amuseDir + ", schedulerType="
                + schedulerType + ", configuration=" + configuration + ", startHub=" + startHub + ", hub=" + hub + "]";
    }

    public Map<String, String> getStatusMap() throws DistributedAmuseException {
        Map<String, String> result = new LinkedHashMap<String, String>();

        result.put("ID", Integer.toString(id));
        result.put("Name", name);
        result.put("Location", location);
        result.put("Gateway", gateway);
        result.put("Amuse dir", amuseDir);
        result.put("Scheduler type", schedulerType);

        result.put("Java path", configuration.getJava());
        result.put("Mpiexec enabled", Boolean.toString(configuration.isMpiexecEnabled()));
        if (configuration.isMpiexecEnabled()) {
            result.put("Mpiexec", configuration.getMpiexec());
        }

        return result;
    }

}
