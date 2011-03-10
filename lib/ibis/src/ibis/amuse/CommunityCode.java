package ibis.amuse;


/**
 * Code representing a code. Also provides the (native) interface to it.
 * 
 * @author Niels Drost
 *
 */
public class CommunityCode {
    
    private final String name;
    
    CommunityCode(String name) {
        this.name = name;
        
        System.loadLibrary("ibis-amuse-" + name);
    }
    
    private native void call() throws Exception;
    
    public Result call(Call call) {
        Result result = new Result(call.getId(), call.getFunctionID(),
                new Exception("cannot actually do any calls yet!"));
        return result;
    }
    
    
    
}
