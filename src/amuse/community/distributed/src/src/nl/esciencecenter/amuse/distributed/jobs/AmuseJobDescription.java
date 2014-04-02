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
package nl.esciencecenter.amuse.distributed.jobs;

import java.io.Serializable;

/**
 * @author Niels Drost
 * 
 */
public abstract class AmuseJobDescription implements Serializable {

    private static final long serialVersionUID = 0L;
    
    private static int nextID = 0;

    private synchronized int getNextID() {
        return nextID++;
    }

    protected final int id;
    protected final String stdoutFile;
    protected final String stderrFile;
    protected final String label;

    public AmuseJobDescription(String stdoutFile, String stderrFile, String label) {
        this.id = getNextID();

        this.stdoutFile = stdoutFile;
        this.stderrFile = stderrFile;
        this.label = label;
    }

    public int getID() {
        return id;
    }

    public String getStdoutFile() {
        return stdoutFile;
    }

    public String getStderrFile() {
        return stderrFile;
    }

    public String getLabel() {
        return label;
    }
    
    public abstract int getNrOfSlots();
    
    public abstract String getType();

}
