import threading
import sys
import logging
import socket


from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface
from amuse.community.interface.common import CommonCode
from amuse.units import units
from amuse.support import options

from distributed_datamodel import Resources, Resource
from distributed_datamodel import Reservations, Reservation

logger = logging.getLogger(__name__)

class OutputHandler(threading.Thread):
    
    def __init__(self, stream, port):
        threading.Thread.__init__(self)
        self.stream = stream

        logging.getLogger("channel").debug("output handler connecting to distributed code")
        
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        
        address = ('localhost', port)
        
        try:
            self.socket.connect(address)
        except:
            raise exceptions.CodeException("Could not connect to distributed code at " + str(address))
        
        self.socket.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
        
        self.socket.sendall('TYPE_OUTPUT'.encode('utf-8'))

        #fetch ID of this connection
        
        result = SocketMessage()
        result.receive(self.socket)
        
        self.id = result.strings[0]
        
        self.daemon = True
        self.start()
        
    def run(self):
        
        while True:
            logging.getLogger("channel").debug("receiving data for output")
            data = self.socket.recv(1024)
            
            if len(data) == 0:
                logging.getLogger("channel").debug("end of output", len(data))
                return
            
            logging.getLogger("channel").debug("got %d bytes", len(data))
            
            self.stream.write(data)

class DistributedAmuseInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
	Distributed Amuse Code
    
        .. [#] The Distributed Amuse project is a collaboration between Sterrewacht Leiden and The Netherlands eScience Center.
    """

    classpath = ['.', 'worker.jar', 'src/dist/*']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="distributed_worker_java", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        
        port = self.get_worker_port()
        
        #logging.basicConfig(level=logging.DEBUG)
        
        logger.debug("running on port %d", port)

#        self.stdoutHandler = OutputHandler(sys.stdout, port)
#        self.stderrHandler = OutputHandler(sys.stderr, port)

        options.GlobalOptions.instance().override_value_for_option("channel_type", "distributed")
        options.GlobalOptions.instance().override_value_for_option("port", port)


    @option(choices=['mpi','remote','distributed', 'sockets'], sections=("channel",))
    def channel_type(self):
        return 'sockets'
    
    @legacy_function
    def get_worker_port():
        """
        Returns the server socket port of the code. Used by the distributed channel
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_webinterface_url():
        """
        Returns the url of the webinterface running inside the distributed code
        """
        function = LegacyFunctionSpecification()
        function.addParameter("address", dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def new_resource():
        """
        Define a new resource. This function returns an index that can be used to refer
        to this resource.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('resource_id', dtype='int32', direction=function.OUT)
        function.addParameter("name", dtype='string', direction=function.IN)
        function.addParameter("location", dtype='string', direction=function.IN)
        function.addParameter("amuse_dir", dtype='string', direction=function.IN)
        function.addParameter("gateway", dtype='string', direction=function.IN, default=[""])
        function.addParameter("scheduler_type", dtype='string', direction=function.IN, default=["ssh"])
        function.addParameter('start_hub', dtype='int32', direction=function.IN, default=1)
        function.addParameter("options", dtype='string', direction=function.IN, default = [""])
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_resource_state():
        """
        Get all the attributes of a resource.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('resource_id', dtype='int32', direction=function.IN)
        function.addParameter("name", dtype='string', direction=function.OUT)
        function.addParameter("location", dtype='string', direction=function.OUT)
        function.addParameter("gateway", dtype='string', direction=function.OUT)
        function.addParameter("amuse_dir", dtype='string', direction=function.OUT)
        function.addParameter("scheduler_type", dtype='string', direction=function.OUT)
        function.addParameter('start_hub', dtype='int32', direction=function.OUT)
        function.addParameter("options", dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_resource():
        """
        Remove the definition of resource from the code. After calling this function the resource is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('resource_id', dtype='int32', direction=function.IN,
            description = "Index of the resource to be removed. This index must have been returned by an earlier call to :meth:`new_resource`")

        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            resource was removed from the model
        -1 - ERROR
            resource could not be removed
        -2 - ERROR
            not yet implemented
        """
        return function
    
    @legacy_function
    def new_reservation():
        """
        Reserve one or more nodes for later use by the simulation.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('reservation_id', dtype='int32', direction=function.OUT)
        function.addParameter("resource_name", dtype='string', direction=function.IN)
        function.addParameter("queue_name", dtype='string', direction=function.IN, default=[""])
        function.addParameter("node_count", dtype='int32', direction=function.IN, default = 1)
        function.addParameter("time", dtype='int32', direction=function.IN, unit = units.minute, default = 60)
        function.addParameter("slots_per_node", dtype='int32', direction=function.IN, default = 1)
        function.addParameter("node_label", dtype='string', direction=function.IN, default = ["default"])
        function.addParameter("options", dtype='string', direction=function.IN, default = [""])
        function.addParameter('count', dtype='int32', direction=function.LENGTH)

        function.result_type = 'int32'
        return function
    
 

    @legacy_function
    def get_reservation_state():
        """
        Get all attributes of a reservation
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('reservation_id', dtype='int32', direction=function.IN)
        function.addParameter("resource_name", dtype='string', direction=function.OUT)
        function.addParameter("queue_name", dtype='string', direction=function.OUT)
        function.addParameter("node_count", dtype='int32', direction=function.OUT)
        function.addParameter("time", dtype='int32', direction=function.OUT, unit = units.minute)
        function.addParameter("slots_per_node", dtype='int32', direction=function.OUT)
        function.addParameter("node_label", dtype='string', direction=function.OUT)
        function.addParameter("options", dtype='string', direction=function.OUT)
        function.addParameter('status', dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)

        function.result_type = 'int32'
        return function

    
    @legacy_function
    def get_reservation_status():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('resource_id', dtype='int32', direction=function.IN)
        function.addParameter('status', dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_reservation():
        """
        Delete a reservation.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('reservation_id', dtype='int32', direction=function.IN)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def wait_for_reservations():
        """
        Wait until all reservations are started, and all nodes are available to run jobs and/or workers
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def submit_script_job():
        """
        Submit a job, specified by a script
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.OUT)
        function.addParameter('script_name', dtype='string', direction=function.IN)
        function.addParameter('arguments', dtype='string', direction=function.IN)
        function.addParameter('script_dir', dtype='string', direction=function.IN)
        function.addParameter("node_label", dtype='string', direction=function.IN, default = ["default"])
        function.addParameter("re_use_code_files", dtype='int32', direction=function.IN, default = 0)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_script_job_state():
        """
        Get all attributes of a script job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter('script_name', dtype='string', direction=function.OUT)
        function.addParameter('arguments', dtype='string', direction=function.OUT)
        function.addParameter('script_dir', dtype='string', direction=function.OUT)
        function.addParameter("node_label", dtype='string', direction=function.OUT)
        function.addParameter("re_use_code_files", dtype='int32', direction=function.OUT)
        function.addParameter("status", dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_script_job_status():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter('status', dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_script_job():
        """
        Delete (cancel) a script job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def wait_for_script_jobs():
        """
        Wait until all script jobs are done.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def submit_function_job():
        """
        Submit a job, specified by a pickle of the function, and a pickle of the arguments.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.OUT)
        function.addParameter('function', dtype='string', direction=function.IN)
        function.addParameter('arguments', dtype='string', direction=function.IN)
        function.addParameter("node_label", dtype='string', direction=function.IN, default = ["default"])
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_function_job_state():
        """
        Get all attributes of a pickled job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter("node_label", dtype='string', direction=function.OUT)
        function.addParameter("status", dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    
    
    @legacy_function
    def get_function_job_status():
        """
        Get all attributes of a pickled job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter("status", dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_function_job_result():
        """
        Get a result of a picked function job. Will block until the result is available
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter('result', dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_function_job():
        """
        Delete (cancel) a script job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_worker_state():
        """
        Get all attributes of a pickled job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('worker_id', dtype='int32', direction=function.IN)

        function.addParameter('executable', dtype='string', direction=function.OUT)
        function.addParameter("node_label", dtype='string', direction=function.OUT)
        function.addParameter("worker_count", dtype='int32', direction=function.OUT)
        function.addParameter("node_count", dtype='int32', direction=function.OUT)
        function.addParameter("thread_count", dtype='int32', direction=function.OUT)
        function.addParameter("status", dtype='string', direction=function.OUT)
        
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_worker_status():
        """
        Get all attributes of a pickled job
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('worker_id', dtype='int32', direction=function.IN)
        function.addParameter("status", dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_number_of_workers():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_workers', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_worker_ids():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index', dtype='int32', direction=function.IN) # probably unused, but required to get 'count'
        function.addParameter('id_of_the_worker', dtype='int32', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_current_error():
        """When a function returns an error, this will retrieve
        a description (if possible)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('string', 
            dtype='string',
            direction=function.OUT,
            description = "description of the error"
        )
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def end_all():
        """
        Stop all jobs, resources and reservations
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    def cleanup_code(self):
        del options.GlobalOptions.instance().overriden_options["channel_type"]
        self.end_all()
        return 0
    
    def new_worker(self):
        raise exceptions.AmuseException("Can't add to 'workers' directly. Create community code instances in the usual way instead.")
    
    def delete_worker(self):
        raise exceptions.AmuseException("Can't remove from 'workers' directly. Stop community code instances in the usual way instead.")
    
    

class DistributedAmuse(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  DistributedAmuseInterface(**options), **options)
    
    def submit_function_job(self, function, *args, **kwargs):
        # pickle the input function
        return self.overridden().submit_function_job(*args, **kwargs)
    
    def get_function_job_result(self, function, *args, **kwargs):
        result = self.overridden().get_function_job_result(*args, **kwargs)
        #un-pickle
        return result
    
    def define_particle_sets(self, object):
        object.define_super_set('items', ['resources', 'reservations', 'script_jobs', 'function_jobs', '_workers'])
        
        #resources
        object.define_set('resources', 'resource_id')
        object.set_new('resources', 'new_resource')
        object.set_delete('resources', 'delete_resource')
        object.add_getter('resources', 'get_resource_state')
        object.mapping_from_name_to_set_definition['resources'].particles_factory = Resources
        
        #reservations
        object.define_set('reservations', 'reservation_id')
        object.set_new('reservations', 'new_reservation')
        object.set_delete('reservations', 'delete_reservation')
        object.add_getter('reservations', 'get_reservation_state')
        object.add_getter('reservations', 'get_reservation_status', names = ('status',))
        object.mapping_from_name_to_set_definition['reservations'].particles_factory = Reservations
        
        #script jobs
        object.define_set('script_jobs', 'job_id')
        object.set_new('script_jobs', 'submit_script_job')
        object.set_delete('script_jobs', 'delete_script_job')
        object.add_getter('script_jobs', 'get_script_job_state')
        object.add_getter('script_jobs', 'get_script_job_status', names = ('status',))
        
        #function jobs
        object.define_set('function_jobs', 'job_id')
        object.set_new('function_jobs', 'submit_function_job')
        object.set_delete('function_jobs', 'delete_function_job')
        object.add_getter('function_jobs', 'get_function_job_state')
        object.add_getter('function_jobs', 'get_function_job_status')
        
        #workers
        object.define_set('_workers', 'worker_id')
        object.set_new('_workers', 'new_worker')
        object.set_delete('_workers', 'delete_worker')
        object.add_getter('_workers', 'get_worker_state')
        object.add_getter('_workers', 'get_worker_status', names = ('status',))
        
    @property
    def workers(self):
        self.update_workers_particle_set()
        return self._workers
    
    def update_workers_particle_set(self):
        """
        Update the "workers" particle set after new instances of codes have been
        created or previously created instances have been stopped.
        """
        old_ids = set(self._workers.get_all_indices_in_store())
        number_of_workers = self.get_number_of_workers()
        if not number_of_workers == 0:
            new_ids = set(self.get_worker_ids(range(number_of_workers)))
        
        ids_to_remove = old_ids - new_ids
        ids_to_add = new_ids - old_ids
        if not len(ids_to_remove) == 0:
            self._workers._remove_indices_in_attribute_storage(list(ids_to_remove))
        if not len(ids_to_add) == 0:
            self._workers._add_indices_in_attribute_storage(list(ids_to_add))
    
