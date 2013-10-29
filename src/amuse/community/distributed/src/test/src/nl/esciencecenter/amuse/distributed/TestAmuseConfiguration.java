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
package nl.esciencecenter.amuse.distributed;

import static org.junit.Assert.*;

import java.io.File;
import org.junit.Test;

/**
 * @author Niels Drost
 * 
 */
public class TestAmuseConfiguration {

    @Test
    public void testAmuseConfiguration() throws Exception {
        File testFile = new File("/home/niels/workspace/amuse");
        AmuseConfiguration config = new AmuseConfiguration(testFile);

        assertEquals(config.getConfigOption("JAVA_ENABLED"), "yes");
    }

    @Test(expected = DistributedAmuseException.class)
    public void testAmuseConfiguration_InvalidOption_DistributedAmuseException() throws Exception {
        File testFile = new File("/home/niels/workspace/amuse/config.mk");
        AmuseConfiguration config = new AmuseConfiguration(testFile);

        assertEquals(config.getConfigOption("SOME_INVALID_CONFIG_OPTION"), "yes");
    }

}
