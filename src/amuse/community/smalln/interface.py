from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

# *** This script, together with the defaults in
# *** GravitationalDynamicsInterface, will be used to generate both
# *** the header file interface.h and the stub interface.cc.

class SmallNInterface(CodeInterface,
                      GravitationalDynamicsInterface):
    """
    Self-contained few-body integrator, using a fourth-order,
    shared-timestep, time-symmetrized Hermite scheme including a
    modified unperturbed treatment of close approaches.
    """

    # Interface specification.

    include_headers = ['interface.h']

    def __init__(self, **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker='smallN_worker',
            **options
        )

    # Interface functions:

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
        function.addParameter('radius', dtype='float64', direction=function.IN,
                 description = "The radius of the particle", default = 0)
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

    # Inheritance from GravitationalDynamicsInterface means that
    # functions in the standard interface don't need to be defined.

    # Additional functions defined here will be reflected in
    # interface.h and must be provided in interface.cc in order for
    # smallN_worker to build.
    @legacy_function
    def set_time():
        """
        Set the model time. Should use set_begin_time
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "The model time to start at", unit = nbody_system.time)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Time value was changed
        -2 - ERROR
            The code does not support setting the  time
        """
        return function
   
        
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
    def set_gamma():
        """
        Set the unperturbed threshold.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gamma():
        """
        Get the unperturbed threshold.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_allow_full_unperturbed():
        """
        Set flag to allow full unperturbed motion.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('allow_full_unpert', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_allow_full_unperturbed():
        """
        Get flag to allow full unperturbed motion.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('allow_full_unpert', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_cm_index():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('cm_index', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_cm_index():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('cm_index', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_break_scale():
        """
        Set the scale at which smallN should stop.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('break_scale', dtype='float64',
                              direction=function.IN, unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_structure_check_interval():
        """
        Set the time scale at which smallN should check structure.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('check_interval', dtype='float64',
                              direction=function.IN, unit = nbody_system.time)
        function.result_type = 'int32'
        return function

    @legacy_function
    def is_over():
        """
        Return 1 iff the run is over, according to analyze().
        """
        function = LegacyFunctionSpecification()
        function.addParameter('over', dtype='int32',
                              direction=function.OUT)
        function.addParameter('rlimit', dtype='float64',
                              direction=function.IN, unit = nbody_system.length,
                              default = 0 )
        function.addParameter('verbose', dtype='int32',
                              direction=function.IN, default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def update_particle_tree():
        """
        Update the internal particle tree to reflect the binary
        structure created in evolve
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_particles_added():
        """
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('n_added', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_id_of_added_particle():
        """
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_update', dtype='int32',
                              direction=function.IN, 
                 description = 'index in the updated particles list')
        function.addParameter('index_of_particle', dtype='int32',
                              direction=function.OUT)
        function.can_handle_array = True
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_children_of_particle():
        """
        Return the number of particles added or deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32',
                              direction=function.IN, 
                 description = 'index of the parent particle',
                 unit = INDEX)
        function.addParameter('child1', dtype='int32', direction=function.OUT,
                description = 'index of the first child particle, -1 if none',
                unit = LINK('particles') )
        function.addParameter('child2', dtype='int32', direction=function.OUT,
                unit = LINK('particles'))
        function.can_handle_array = True
        function.result_type = 'int32'
        return function

class SmallN(GravitationalDynamics):

    # The actual module.

    def __init__(self, convert_nbody = None, **keyword_arguments):
        legacy_interface = SmallNInterface(**keyword_arguments)

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
        #     SmallN.parameters.timestep_parameter = xxx

        object.add_method_parameter(
            "get_eta",                   # getter name in interface.cc
            "set_eta",                   # setter name in interface.cc
            "timestep_parameter",        # python parameter name
            "timestep parameter",        # description
            default_value = 0.14
        )
        
        object.add_method_parameter(
            "get_gamma",                 # getter name in interface.cc
            "set_gamma",                 # setter name in interface.cc
            "unperturbed_threshold",     # python parameter name
            "unperturbed threshold",     # description
            default_value = 1.e-6
        )
        
        object.add_method_parameter(
            "get_allow_full_unperturbed",  # getter name in interface.cc
            "set_allow_full_unperturbed",  # setter name in interface.cc
            "allow_full_unperturbed",      # python parameter name
            "full unperturbed motion",     # description
            default_value = 1
        )
        
        object.add_method_parameter(
            "get_cm_index",		   # getter name in interface.cc
            "set_cm_index",		   # setter name in interface.cc
            "cm_index",			   # python parameter name
            "current CM index",	           # description
            default_value = 100000
        )
        
        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )
        
        
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        object.add_getter("particles", 'get_children_of_particle')

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        # Turn interface functions into methods.

        object.add_method("new_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                object.NO_UNIT
            ),
            (
                object.INDEX,
                object.ERROR_CODE
            )
        )

    def update_particle_set(self):
        """
        update the particle set after the code has added binaries
              
        """
        
        number_of_updated_particles = self.get_number_of_particles_added()
        #print "number_of_updated_particles =", number_of_updated_particles
        
        if number_of_updated_particles == 0:
            return
        
        indices_in_update_list = range(number_of_updated_particles)
        indices_to_add = self.get_id_of_added_particle(indices_in_update_list)
        #print "indices_to_add:", indices_to_add, indices_in_update_list
        
        incode_storage = self.particles._private.attribute_storage
        
        if len(indices_to_add) > 0:
            incode_storage._add_indices(indices_to_add)
