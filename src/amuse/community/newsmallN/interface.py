from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

# *** This script, together with the defaults in
# *** GravitationalDynamicsInterface, will be used to generate both
# *** the header file interface.h and the stub interface.cc.

class smallNInterface(CodeInterface,
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

    # Inheritance from GravitationalDynamicsInterface means that
    # functions in the standard interface don't need to be defined.

    # Additional functions defined here will be reflected in
    # interface.h and must be provided in interface.cc in order for
    # smallN_worker to build.

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
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='float64',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gamma():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_allow_full_unperturbed():
        """
        Set the current time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('allow_full_unpert', dtype='int32',
                              direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_allow_full_unperturbed():
        """
        Set the current system time step parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('allow_full_unpert', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function

class smallN(GravitationalDynamics):

    # The actual module.

    def __init__(self, convert_nbody = None, **keyword_arguments):
        legacy_interface = smallNInterface(**keyword_arguments)

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
        #     smallN.parameters.timestep_parameter = xxx

        object.add_method_parameter(
            "get_eta",                   # getter name in interface.cc
            "set_eta",                   # setter name in interface.cc
            "timestep_parameter",        # python parameter name
            "timestep parameter",        # description
            default_value = 0.14 | units.none
        )
        
        object.add_method_parameter(
            "get_gamma",                 # getter name in interface.cc
            "set_gamma",                 # setter name in interface.cc
            "unperturbed_threshold",     # python parameter name
            "unperturbed threshold",     # description
            default_value = 1.e-6 | units.none
        )
        
        object.add_method_parameter(
            "get_allow_full_unperturbed",  # getter name in interface.cc
            "set_allow_full_unperturbed",  # setter name in interface.cc
            "allow_full_unperturbed",      # python parameter name
            "full unperturbed motion",     # description
            default_value = 1 | units.none
        )
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        # Turn interface functions into methods.

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
                units.none
            ),
            (
                object.INDEX,
                object.ERROR_CODE
            )
        )

        object.add_method("get_eta", (), (units.none, object.ERROR_CODE))
        object.add_method("set_eta", (units.none), (object.ERROR_CODE))
        object.add_method("get_gamma", (), (units.none, object.ERROR_CODE))
        object.add_method("set_gamma", (units.none), (object.ERROR_CODE))
        object.add_method("get_allow_full_unperturbed",
                          (), (units.none, object.ERROR_CODE))
        object.add_method("set_allow_full_unperturbed",
                          (units.none), (object.ERROR_CODE))
