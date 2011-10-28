package ibis.amuse;

import java.util.Map;

import ibis.ipl.IbisIdentifier;

public interface MPIProfilingCollectorMBean {
	
	public Object getSentBytesPerIbis();
	
	public Object getReceivedBytesPerIbis();
}
