from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class FDPSInterface(
        CodeInterface,
        LiteratureReferencesMixIn,
        GravitationalDynamicsInterface,
        #StoppingConditionInterface,
        #SinglePointGravityFieldInterface,
        ):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(
                self,
                name_of_the_worker="fdps_sd_worker", 
                **options)
    
#    @legacy_function
#    def get_number_of_particles():
#        function = LegacyFunctionSpecification()  
#        function.can_handle_array = True 
#        function.addParameter('value', dtype='int32', direction=function.OUT)
#        function.result_type = 'int32'
#        return function    
    @legacy_function
    def set_time_step():
        """
        Update timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.IN,
            description = "timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_epsilon_squared():
        """
        Get epsilon^2, a softening parameter for gravitational potentials with point particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.OUT,
            description = "epsilon^2, a softening parameter for gravitational potentials with point particles",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_epsilon_squared():
        """
        Set epsilon^2, a softening parameter for gravitational potentials with point particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.IN,
            description = "epsilon^2, a softening parameter for gravitational potentials with point particles",
            unit = nbody_system.length * nbody_system.length)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_theta_for_tree():
        """
        Get theta, the opening angle for building the tree: between 0 and 1.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('theta_for_tree', dtype='float64', direction=function.OUT,
            description = "theta, the opening angle for building the tree: between 0 and 1")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_theta_for_tree():
        """
        Set theta, the opening angle for building the tree: between 0 and 1.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('theta_for_tree', dtype='float64', direction=function.IN,
            description = "theta, the opening angle for building the tree: between 0 and 1")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_group_limit_for_tree():
        """
        Get group limit.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('group_limit', dtype='int32', direction=function.OUT,
            description = "group limit")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_group_limit_for_tree():
        """
        Set group limit.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('group_limit', dtype='int32', direction=function.IN,
            description = "group limit")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_leaf_limit_for_tree():
        """
        Get leaf limit.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('leaf_limit', dtype='int32', direction=function.OUT,
            description = "leaf limit")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_leaf_limit_for_tree():
        """
        Set leaf limit.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('leaf_limit', dtype='int32', direction=function.IN,
            description = "leaf limit")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    
class FDPS(GravitationalDynamics, GravityFieldCode):

    __interface__ = FDPSInterface

    def __init__(self, convert_nbody = None, **options):

        legacy_interface = self.__interface__(**options)

        GravitationalDynamics.__init__(
                self,
                legacy_interface,
                convert_nbody,
                **options
                )

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        GravityFieldCode.define_state(self, object)
        #self.stopping_conditions.define_state(object)

    def define_parameters(self, object):
        object.add_method_parameter(
                "get_begin_time",
                "set_begin_time",
                "begin_time",
                "model time to start the simulation at",
                default_value = 0.0 | nbody_system.time
                )

        object.add_method_parameter(
            "get_epsilon_squared",
            "set_epsilon_squared",
            "epsilon_squared",
            "smoothing parameter for gravity calculations", 
            default_value = 0.125 | nbody_system.length * nbody_system.length
        )
        
        object.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "opening_angle", 
            "opening angle, theta, for building the tree: between 0 and 1", 
            default_value = 0.75
        )

        object.add_method_parameter(
            "get_group_limit_for_tree",
            "set_group_limit_for_tree",
            "group_limit", 
            "group limit",
            default_value = 64
        )

        object.add_method_parameter(
            "get_leaf_limit_for_tree",
            "set_leaf_limit_for_tree",
            "leaf_limit", 
            "leaf limit",
            default_value = 8
        )
        
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            default_value = (1.0/64.0)| nbody_system.time
        )        

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            "get_time",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_time_step",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_time_step",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )

        
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        
        #self.stopping_conditions.define_particle_set(object)
