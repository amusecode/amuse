from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

# *** This script, together with the defaults in
# *** GravitationalDynamicsInterface, will be used to generate both
# *** the header file interface.h and the stub interface.cc.  Since
# *** interface.cc has been  hand-coded to implement the details,
# *** MAKE SURE TO SAVE IT SOMEWHERE, as build.py can overwrite it!

class ph4Interface(CodeInterface,
                   GravitationalDynamicsInterface):
    """
    Parallel, GPU-accelerated, N-body integration module with block
    time steps, using a 4th-order Hermite integration scheme.
    """

    # Interface specification.

    include_headers = ['interface.h']

    MODE_GPU = 'gpu'
    MODE_CPU = 'cpu'
    
    def __init__(self, mode = MODE_CPU, **options):
        print "worker code =", self.name_of_the_muse_worker(mode)
        CodeInterface.__init__(
            self,
            name_of_the_worker=self.name_of_the_muse_worker(mode),
            **options
        )
    
    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_CPU:
            return 'ph4_worker'
        elif mode == self.MODE_GPU:
            return 'ph4_worker_gpu'
        else:
            return 'ph4_worker'
        
    # Inheritance from GravitationalDynamicsInterface means that
    # functions in the standard interface don't need to be defined.
    # See interface.py.2 for a laboriously hand-coded version written
    # before I discovered this fact!    (Steve McMillan, 10/10)

    # Additional functions defined here will be reflected in
    # interface.h and must be provided in interface.cc in order for
    # ph4_worker to build.

    # The following functions aren't defined in the default interface:

    @legacy_function
    def set_eta():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_gpu():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('gpu', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gpu():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('gpu', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time():
        """
        Set the current system time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('system_time', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_number_of_particles_updated():
        """
        Return the number of particles added or deleted during the last evolve
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_id_of_updated_particle():
        """
        Return the number of particles added or deleted during the last evolve
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_update', dtype='int32',
                              direction=function.IN, 
                              description = 'index in the updated particles list')
        function.addParameter('index_of_particle', dtype='int32',
                              direction=function.OUT)
        function.addParameter('kind_of_update', dtype='int32',
                              direction=function.OUT,
                              description = 'kind of update (2, addition), (1, deletion)')
        function.can_handle_array = True
        function.result_type = 'int32'
        return function


class ph4(GravitationalDynamics):

    # The actual module.

    def __init__(self, convert_nbody = None, **keyword_arguments):
        legacy_interface = ph4Interface(**keyword_arguments)

        GravitationalDynamics.__init__(self,
                                       legacy_interface,
                                       convert_nbody,
                                       **keyword_arguments)

    def define_parameters(self, object):

        # Set/get parameters specific to the module, not part of the
        # standard interface.  Accessors used here must be defined
        # above and reflected in interface.cc.  Python access is
        # (e.g.)
        #
        #	ph4.parameters.timestep_parameter = xxx

        object.add_method_parameter(
            "get_eta",			# getter name in interface.cc
            "set_eta",			# setter name in interface.cc
            "timestep_parameter",	# python parameter name
            "timestep parameter",	# description
            units.none,			# units
            0.14 | units.none		# default
        )

        object.add_method_parameter(
            "get_eps2",			# already defined in standard interface
            "set_eps2",			# already defined in standard interface
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.01 | nbody_system.length * nbody_system.length
        )

        object.add_method_parameter(
            "get_gpu",			# getter name in interface.cc
            "set_gpu",			# setter name in interface.cc
            "use_gpu",			# python parameter name
            "use GPU",			# description
            units.none,			# units
            1 | units.none		# default
        )

    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.support.data.incode_storage module, as
        it uses a lot of internal methods and info!
        
        """
        number_of_updated_particles, error \
		= self.get_number_of_particles_updated()
        
        if number_of_updated_particles == 0:
            return
        
        indices_in_update_list = range(number_of_updated_particles)
        particle_indices, updates, erros \
		= self.get_id_of_updated_particle(indices_in_update_list)
        
        incode_storage = self.particles._private.attribute_storage
        
        indices_to_remove = []
        indices_to_add = []
        for index, status in zip(particle_indices, updates):
            if status == 1:			# deletion
                indices_to_remove.append(index)
            elif status == 2:			# addition
                indices_to_add.append(index)

        print "indices_to_remove:", indices_to_remove
        print "indices_to_add:", indices_to_add
        
        if len(indices_to_remove) > 0:
            incode_storage._remove_indices(indices_to_remove)
        
        if len(indices_to_add) > 0:
            incode_storage._add_indices(indices_to_add)
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        # Similarly, we can add module-specific methods, if desired.
        # See hermite0/interface.py for examples.

