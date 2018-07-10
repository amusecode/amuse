==================
Distributed AMUSE
==================

It is possible to run AMUSE on multiple machines simultaneously. 
The AMUSE script itself always runs on a users' local machine, while workers for codes can be "send out" to remote machines such as workstations, clusters, etc.


Installation
------------

Deploying workers one remote machines requires a full installation of AMUSE on each machine. For each code "sockets" support needs to be present. This is done by default, and should be available for all codes. Note that older versions of Fortran lack the necessary features to support sockets workers.

On each machine, the distributed code also needs to be build. Distributed AMUSE requires a Java Development Kit (JDK), preferably Oracle Java version 7 or 8. The ``configure`` script tries to locate the JDK, but you may need to specify it by hand. For details, see:

.. code-block:: sh

	> ./configure --help

To build distributed amuse run the following at the amuse root:
	
.. code-block:: sh

	> make distributed.code

To check if the installation is set-up properly, run all the tests related to the worker interface:

.. code-block:: sh

	> cd $AMUSE_DIR
	> nosetests -v test/codes_tests/test*implementation.py
	
Note that Distributed AMUSE is mostly tested with the version of MPI includes in the amuse "prerequisites". 
If you get MPI errors while running remote (parallel) workers, try using the install.py script included in AMUSE to install the prerequisites.  

Overview
--------

Usage of Distributed Amuse is (by design) very close to the usage of any other code in AMUSE. 
The main difference being it contains resources, pilots, and jobs, instead of particles.

Resource
	Description of a machine capable of running jobs. 
    For each resource distributed AMUSE will launch a support process (HUB) to facilitate communication and coordination between the workers and AMUSE
	
Pilot
	Process running on a resource (often within a reservation on a resource) waiting for jobs to run on behalf of AMUSE. 
    Can consist of multiple machines.
	
Job
	Worker process.
    Will search for a suitable pilot to run on.

In general, a user will first define resources, then deploy pilots on these resources, and finally create codes that make use of the machines offered by the pilots.


Initializing the Distributed AMUSE code
---------------------------------------

Distributed Amuse can be initialized like any other code:

    >>> from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
    >>> from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
    >>> 
    >>> #initialize code, print output of code to console
    >>> instance = DistributedAmuse(redirection='none')


Parameters
----------

Distributed AMUSE supports a few parameters to adjust settings. 
All parameters need to be set before any resource, pilot or job is made to have effect.

Overview of settings:

debug
	Boolean parameters, defaults to False. 
    If true/enabled, will output additional debugging information and logs, both in the code output, and in a `distributed-amuse-logs` folder on any target machine used.
webinterface_port
	Port on which a simple webinterface is available for monitoring. 
    Defaults to "0", for a port determined automatically.
start_hubs
	To facilitate communication across different networks (with for instance firewalls), as hub is by default started on each resource. 
    This can be turned off if needed, for instance if all resources are within the same network.
worker_queue_timeout
	The user is responsible for making sure enough slots are available to run a worker. 
    If not, it will end up in the queue. 
    The time the worker will wait before giving up can be set using this parameter.
worker_startup_timeout
	The distributed code starts AMUSE workers running the actual codes. 
    This can take a while on some machines. 
    If needed, this parameter can be used to increase the time waited.
	
    >>> instance.parameters.debug = True
    >>> instance.parameters.webinterface_port = 5555
    >>> instance.commit_parameters()
    >>>
    >>> print instance.parameters.webinterface_port
    

Monitoring
----------

Distributed Amuse has a small build-in webinterface for monitoring. 
A utility function is available to get the url:

    >>> import webbrowser
    >>>
    >>> webbrowser.open(instance.get_webinterface_url())

Specifying resources
--------------------

In order to use a remote machine, AMUSE needs to have some information about this resource such as the host name, type of machine, username to gain access, etc.
This can be specified by creating a "Resource" in Distributed AMUSE. 
As a side effect, a communication hub is also started on the (frontend of) the resource.

    >>> resource = Resource()
    >>> resource.name = "some.resource"
    >>> resource.location = "user@machine.example.com"
    >>> resource.scheduler = "ssh"
    >>> resource.amuse_dir = "/home/user/amuse"
    >>>
    >>> instance.resources.add_resource(resource)

Overview of all options:

name
	Some user chosen name for the resource
location
	Address of the resource. Usually a hostname (e.g. somehost.somedomain.com). Could also be an IP address
amuse_dir
	Location of amuse on the remote machine (e.g. /home/user/amuse-svn)
tmp_dir
	Where all temporary files will be put on the remote machine
gateway
	Sometimes a machine is not reachable directly due to firewalls and such. Use this setting to provide an intermediate resource to route traffic via. This resource should already have been created.
scheduler_type
	The type of scheduler present on the remote machine. Defaults to 'ssh' useful for single machines. Current supported scheduler types: 'ssh', 'sge', 'slurm'
hub_queue_name
	Normally the support process is started on the front end. However, it can also be submitted to a queue by specifying it here.
hub_time_minutes
	When a hub is submitted, this option denotes the time the hub will be available.


Starting Pilots
---------------

The next step in running jobs remotely is to start a so-called pilot job on the resource specified previously. This pilot will submit a job to the resource, create necessary communication channels with the main amuse application, and wait for jobs to be started (currently mostly workers)

Note that pilots may not be started for a while. A function is available to wait until all created pilots have started.

    >>> pilot = Pilot()
    >>> pilot.resource_name='local'
    >>> pilot.node_count=1
    >>> pilot.time= 2|units.hour
    >>> pilot.slots_per_node=22
    >>> pilot.label='local'
    >>>
    >>> instance.pilots.add_pilot(pilot)
    >>> 
    >>> print "Pilots:"
    >>> print instance.pilots
    >>> 
    >>> print "Waiting for pilots"
    >>> instance.wait_for_pilots()

Overview of all options:

resource_name
	name of the resource to start the pilot on
queue_name
	queue to use to run the pilot (cluster specific, not used in case of ssh)
node_count
	number of nodes to start the pilot on
time
	time to keep the pilot active
slots_per_node
	number of workers to start on a node. Usually the number of cores, but could be less if memory is a limiting factor, or workers are multi-core capable
label
	label to attach to the pilot. Can be used when starting workers to run workers on specific pilots
options
	Additional options. Usually not required.


Starting jobs
-------------

When running remote workers, they can be started as normal. 
However, AMUSE needs to be signalled to use the distributed code to start them instead of the normal process.
A function is available to enable and disable this.

    >>> print "starting all workers using the distributed code"
    >>> instance.use_for_all_workers()

    >>> print "not using distributed workers any longer"
    >>> instance.use_for_all_workers(enable=False)

Alternatively, you can also explicitly enable the distributed code per worker

    >>> print "using this distributed instance for all distributed workers"
    >>> instance.use_for_all_distributed_workers(enable=True)
    >>> worker = Hermite(channel_type='distributed')

Or, even pass the instance of the distributed code you would like to use, in the rare case you have multiple distributed codes

    >>> worker = Hermite(channel_type='distributed', distributed_instance=instance)

Worker options
--------------

This section lists all the relevant worker options for Distributed AMUSE. 
Most are new, some are also supported in the other channel implementations.
You are normally not required to use any options.

number_of_workers
	Number of worker processes started (thus working as normally the case). 
    Each worker takes up a slot of the pilot (see above)
label
	Label of the pilot to use. By default any pilot with enough free slots found will be used to start this worker. 
    Using the labels an explicit selection can be done.
number_of_threads
	Number of threads used in the process. 
    This can be used to explicitly set the OMP_NUM_THREADS environment variable in the worker
channel_type
	Set this to "distributed" to start workers using the distributed code. 
    Alternatively, use the use_for_all_workers functions as described above to set this by default
distributed_instance
	This is a reference to the distributed instance used to start the worker, in the rare case you have multiple distributed codes.
dynamic_python_code
	Boolean option stating if this code is a dynamic python code. 
    If so, all .py files in the worker directory will be copied to the remote machine before starting the code.


Labels
------

By default workers are started on any available pilot with enough slots available. 
However, sometimes you would like to have more control over which worker is started where, for instance if special hardware is present on some machines.

The concept of labels can be used within Distributed AMUSE to get this functionality.
If a label is attached to a worker (one of the parameters when starting a worker, see above), only pilots with exactly the same label (specified when the pilot is started) are considered candidates for running the worker. 
The name of labels is completely up to the user.

For instance, say a simulation uses a number of workers running on a CPU, and a single GPU worker.
The following code will put all the cpu workers on one machine, and the single gpu worker on another.

    >>> cpu_pilot = Pilot()
    >>> cpu_pilot.resource_name='machine1'
    >>> cpu_pilot.node_count=1
    >>> cpu_pilot.time= 2|units.hour
    >>> cpu_pilot.slots_per_node=30
    >>> cpu_pilot.label='CPU'
    >>> instance.pilots.add_pilot(cpu_pilot)
    >>>
    >>> gpu_pilot = Pilot()
    >>> gpu_pilot.resource_name='machine2'
    >>> gpu_pilot.node_count=1
    >>> gpu_pilot.time= 2|units.hour
    >>> gpu_pilot.slots_per_node=1
    >>> gpu_pilot.label='GPU'
    >>> instance.pilots.add_pilot(gpu_pilot)
    >>>
    >>> ...
    >>> worker1 = Hermite(label='CPU')
    >>> worker2 = Bonsai(label='GPU')
    >>>
    >>> #will not start due to a lack of slots.
    >>> worker3 = Bonsai(label='GPU')
 

Examples
--------

AMUSE contains a number of examples for the distributed code. See examples/applications/

Gateways
--------

Gateways can be used in case of connectivity problems between machines, such as firewalls and private IP addresses. 
This is for instance the case at the LGM. 
A gateway is started like any other resource (and thus require a valid installation of AMUSE on each gateway). 
This resource can then be specified to be a "gateway" to another resource. 
In this case all ssh connections will be made via the gateway, so make sure you can login from the gateway to the target machine without using a password, as well as from your local machine.

Commonly Encountered Problems
-----------------------------

Most initial setup problems with the Distributed AMUSE code can be solved by checking:

- Can you login to each machine you plan to use using ssh without using a password? 
  See for instance here on how to set this up: https://www.thegeekstuff.com/2008/11/3-steps-to-perform-ssh-login-without-password-using-ssh-keygen-ssh-copy-id/
- Did you configure a Java JDK version 1.7 or higher using ./configure? 
  Check the content of config.mk to see which java is used, and what version was detected. 
  Make sure to do a "make clean" and "make" in case you make any changes. This should also be done on all machines.
- Is AMUSE configured properly on each and every machine? 
  Running the code implementation tests is a good way of spotting issues:

    >>> nosetests -v test/codes_tests/test_*_implementation.py

- Are the settings provided for each resource correct (username, amuse location, etc)
- Have you set the correct mpiexec in ./configure? This setting is normally not used by AMUSE, so you may only now notice it is misconfigured

In case this does not help, it is probably best to check the output for any errors. 
Normally worker output is discarded by most scripts. 
Use 'redirect=none' to see the output of the workers, a lot of errors show up in this output only. 
There is also a "debug" parameter in Distributed Amuse.
If enabled, output for each pilot will be in a "distributed-amuse-logs" folder in the home of each remote machine used, and additional information is printed to the log from the local AMUSE script.

