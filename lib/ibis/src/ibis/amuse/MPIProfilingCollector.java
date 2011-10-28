package ibis.amuse;

import ibis.ipl.IbisIdentifier;

import java.lang.management.ManagementFactory;
import java.util.HashMap;
import java.util.Map;

import javax.management.MBeanServer;
import javax.management.ObjectName;

public class MPIProfilingCollector extends Thread implements MPIProfilingCollectorMBean  {
	
	private final Map<IbisIdentifier, Long> sent;
	private final Map<IbisIdentifier, Long> received;
	
	private final IbisIdentifier[] ibisList;
	
	public MPIProfilingCollector(IbisIdentifier[] ibisList) throws Exception {
		this.ibisList = ibisList;
		
		sent = new HashMap<IbisIdentifier, Long>();
		received = new HashMap<IbisIdentifier, Long>();
		
		MBeanServer mbs = ManagementFactory.getPlatformMBeanServer(); 
	    ObjectName name = new ObjectName("ibis.amuse:type=MPIProfilingCollector"); 
	    mbs.registerMBean(this, name); 
	    
	    this.setDaemon(true);
	    this.setName(this.getClass().getName());
	    this.start();
	}

	@Override
	public synchronized Map<IbisIdentifier, Long> getSentBytesPerIbis() {
		return new HashMap<IbisIdentifier, Long>(sent);
	}

	@Override
	public synchronized Map<IbisIdentifier, Long> getReceivedBytesPerIbis() {
		return new HashMap<IbisIdentifier, Long>(received);
	}

	public void run() {
		long i = 1000;
		while(true) {
			i += 1000;
			for (IbisIdentifier ibis: ibisList) {
				sent.put(ibis, i);
				received.put(ibis, i);
			}
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				return;
			}
		}
	}
	

}
