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
        
    @legacy_function
    def new_particle():
        """
        Define a new particle in the stellar dynamics code. The
        particle is initialized with the provided mass, radius,
        position and velocity. This function returns an index that
        can be used to refer to this particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32',
                              direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )

        function.addParameter('mass', dtype='float64', direction=function.IN,
                 description = "The mass of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN,
                 description = "The radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN,
                 description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN,
                 description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN,
                 description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN,
                 description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN,
                 description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN,
                 description = "The initial velocity vector of the particle")
        function.addParameter('id', dtype='int32', direction=function.IN,
                 description = "Identifier of the particle, "
                               +"option for restoring state after loading",
                              default = -1)
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            particle was created and added to the model
        -1 - ERROR
            particle could not be created"""
        return function
    
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
    def set_manage_encounters():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('manage_encounters', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_manage_encounters():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('manage_encounters', dtype='int32',
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
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_id_of_updated_particle():
        """
        Return the number of particles added or deleted during the last evolve.
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

    @legacy_function
    def get_binary_energy():
        """
        Return the total energy in all binaries.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'binary_energy', 
            dtype='float64',
            direction=function.OUT
        )
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
        
        object.add_method_parameter(
            "get_manage_encounters",	# getter name in interface.cc
            "set_manage_encounters",	# setter name in interface.cc
            "manage_encounters",	# python parameter name
            "manage close encounters",	# description
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

        print ''
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

        object.add_method("new_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                units.none,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )

        object.add_method("get_binary_energy", (),
	    (nbody_system.energy, object.ERROR_CODE))
