package ibis.amuse;

import java.io.File;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.util.TypedProperties;

public class AmuseConfigOptions {

    private static final Logger logger = LoggerFactory.getLogger(AmuseConfigOptions.class);

    private final File amuseHome;

    private final TypedProperties properties;

    AmuseConfigOptions(String amuseHome) {
        this.amuseHome = new File(amuseHome);
        this.properties = new TypedProperties();

        properties.loadFromFile(amuseHome + "/config.mk");

        logger.info("Amuse home = " + amuseHome);

        if (logger.isDebugEnabled()) {
            for (Object key : properties.keySet()) {
                logger.debug("key: " + key + " value = " + properties.get(key));
            }
        }
    }

    String getPytonPath() {
        return (new File(amuseHome, "test")) + ":" + (new File(amuseHome, "src"));
    }

    File getAmuseHome() {
        return amuseHome;
    }

    public String getMpiexec() {
        return properties.getProperty("MPIEXEC");
    }

    public boolean isMpiexecEnabled() {
        return properties.getBooleanProperty("MPIEXEC_ENABLED");
    }

}
