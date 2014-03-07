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

import nl.esciencecenter.amuse.distributed.jobs.WorkerDescription;

import org.junit.Test;

/**
 * @author Niels Drost
 * 
 */
public class DistributedAmuseIT {

    @Test
    public void test01_Constructor() throws DistributedAmuseException {
        DistributedAmuse da = new DistributedAmuse("/home/niels/workspace/amuse/src/amuse/community/distributed",
                "/home/niels/workspace/amuse", 8678, true);

        da.pilots().newPilot("local", "unlimited", 1, 10, 1, "default", "");
        da.pilots().waitForAllPilots();
        
        WorkerDescription description = new WorkerDescription("some.id", "some.dir/some.worker.executable", "stdout", "stderr","default",
                1,1, 60);
        
        da.jobs().submitWorkerJob(description);
        
        try {
            Thread.sleep(300000);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        da.end();
    }

    @Test
    public void test02_addResources() throws DistributedAmuseException {
        DistributedAmuse da = new DistributedAmuse("/home/niels/workspace/amuse/src/amuse/community/distributed",
                "/home/niels/workspace/amuse", 8678, true);

        da.resources().newResource("some.name", "niels@fs0.das4.cs.vu.nl", null, "/home/niels/amuse", "sge", true, "");
        //da.resourceManager().newResource("lgm", "niels@node04", "niels@fs.lgm.liacs.nl", "/home/niels/amuse", "sge", true);

        da.pilots().newPilot("some.name", "all.q", 1, 10, 1, "default", "");
        
        da.pilots().waitForAllPilots();
        
        try {
            Thread.sleep(60000);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        da.end();
    }

}
