package nl.esciencecenter.amuse.distributed.jobs;

public class FunctionJobDescription extends AmuseJobDescription {
    
    private static final long serialVersionUID = 1L;
    
    private final String function;
    private final String arguments;
    private final String kwarguments;


    public FunctionJobDescription(String function, String arguments, String kwarguments, String stdoutFile, String stderrFile, String nodeLabel) {
        super(stdoutFile, stderrFile, nodeLabel);
        
        this.function = function;
        this.arguments = arguments;
        this.kwarguments = kwarguments;
    }

    @Override
    public int getNrOfSlots() {
        return 1;
    }

    @Override
    public String getType() {
        return "function";
    }
    
    public String getFunction() {
        return function;
    }
    
    public String getArguments() {
        return arguments;
    }
    
    public String getKwArguments() {
        return kwarguments;
    }


    @Override
    public String toString() {
        return "FunctionJobDescription [function length=" + function.length() + ", arguments length=" + arguments.length() + ", kwarguments length=" + kwarguments.length()
                + ", id=" + id + ", stdoutFile=" + stdoutFile + ", stderrFile=" + stderrFile + ", label=" + label + "]";
    }

}
