package ibis.amuse;

import java.io.Serializable;

/**
 * Settings for a worker.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerInfo implements Serializable {

    private static final long serialVersionUID = 1L;

    private String id;

    private int nrOfNodes;

    private int nrOfProcesses;
    
    private int nrOfThreads;

    private String codeName;

    private String codeDir;

    private String amuseHome;

    private String interfaceType;

    private String stdoutFile;

    private String stderrFile;
    
    private boolean isPythonCode;
    
    private boolean copyWorkerCode;

    String getID() {
        return id;
    }

    void setWorkerID(String workerID) {
        this.id = workerID;
    }

    int getNrOfNodes() {
        return nrOfNodes;
    }

    void setNrOfNodes(int nrOfNodes) {
        this.nrOfNodes = nrOfNodes;
    }

    int getNrOfProcesses() {
        return nrOfProcesses;
    }

    void setNrOfProcesses(int nrOfProcesses) {
        this.nrOfProcesses = nrOfProcesses;
    }
    
    int getNrOfThreads() {
        return nrOfThreads;
    }

    void setNrOfThreads(int nrOfThreads) {
        this.nrOfThreads = nrOfThreads;
    }

    String getCodeName() {
        return codeName;
    }

    void setCodeName(String codeName) {
        this.codeName = codeName;
    }

    String getCodeDir() {
        return codeDir;
    }

    void setCodeDir(String codeDir) {
        this.codeDir = codeDir;
    }

    String getAmuseHome() {
        return amuseHome;
    }

    void setAmuseHome(String amuseHome) {
        this.amuseHome = amuseHome;
    }

    String getInterfaceType() {
        return interfaceType;
    }

    void setInterfaceType(String interfaceType) {
        this.interfaceType = interfaceType;
    }

    String getStdoutFile() {
        return stdoutFile;
    }

    void setStdoutFile(String stdoutFile) {
        this.stdoutFile = stdoutFile;
    }

    String getStderrFile() {
        return stderrFile;
    }

    void setStderrFile(String stderrFile) {
        this.stderrFile = stderrFile;
    }

    boolean isPythonCode() {
        return isPythonCode;
    }
    
    void setIsPythonCode(boolean isPythonCode) {
        this.isPythonCode = isPythonCode;
    }
    
    boolean copyWorkerCode() {
        return copyWorkerCode;
    }
    
    void setCopyWorkerCode(boolean copyWorkerCode) {
        this.copyWorkerCode = copyWorkerCode;
    }
}
