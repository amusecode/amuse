from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics

class BHTreeInterface(CodeInterface, LiteratureReferencesMixIn, GravitationalDynamicsInterface, StoppingConditionInterface):
    """
        .. [#] Barnes, J., Hut, P., A Hierarchical O(N log N) force-calculation algorithm, *Nature*, **4**, 324 (1986)   
    """
    include_headers = ['bhtree_code.h', 'worker_code.h', 'stopcond.h']
    
    def __init__(self, convert_nbody = None, **kwargs):
        CodeInterface.__init__(self, name_of_the_worker="bhtree_worker", **kwargs)
        """
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
            
        self.convert_nbody = convert_nbody
        """
        LiteratureReferencesMixIn.__init__(self)

    
       
    
    
    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
    
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
            description = "epsilon^2, a softening parameter for gravitational potentials with point particles")
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
            description = "epsilon^2, a softening parameter for gravitational potentials with point particles")
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
    def get_use_self_gravity():
        """
        Get use_self_gravity flag, the flag for usage of self gravity, 1 or 0 (true or false).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('use_self_gravity', dtype='int32', direction=function.OUT,
            description = "flag for usage of self gravity, 1 or 0 (true or false)")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_use_self_gravity():
        """
        Set use_self_gravity flag, the flag for usage of self gravity, 1 or 0 (true or false).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('use_self_gravity', dtype='int32', direction=function.IN,
            description = "flag for usage of self gravity, 1 or 0 (true or false)")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_ncrit_for_tree():
        """
        Get Ncrit, the maximum number of particles sharing an interaction list.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('ncrit_for_tree', dtype='int32', direction=function.OUT,
            description = "Ncrit, the maximum number of particles sharing an interaction list")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_ncrit_for_tree():
        """
        Set Ncrit, the maximum number of particles sharing an interaction list.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('ncrit_for_tree', dtype='int32', direction=function.IN,
            description = "Ncrit, the maximum number of particles sharing an interaction list")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_dt_dia():
        """
        Get the time interval between diagnostics output.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64', direction=function.OUT,
            description = "the time interval between diagnostics output")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_dt_dia():
        """
        Set the time interval between diagnostics output.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64', direction=function.IN,
            description = "the time interval between diagnostics output")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
class BHTree(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)
        
        legacy_interface = BHTreeInterface(**options)
        
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )     
            
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_epsilon_squared",
            "set_epsilon_squared", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.125 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            default_value = 0.015625 | nbody_system.time
        )
        object.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "opening_angle", 
            "opening angle, theta, for building the tree: between 0 and 1", 
            default_value = 0.75 | units.none
        )
        object.add_method_parameter(
            "get_use_self_gravity",
            "set_use_self_gravity",
            "use_self_gravity", 
            "flag for usage of self gravity, 1 or 0 (true or false)", 
            default_value = 1 | units.none
        )
        object.add_method_parameter(
            "get_ncrit_for_tree",
            "set_ncrit_for_tree",
            "ncrit_for_tree", 
            "Ncrit, the maximum number of particles sharing an interaction list", 
            default_value = 12 | units.none
        )
        object.add_method_parameter(
            "get_dt_dia",
            "set_dt_dia",
            "dt_dia", 
            "time interval between diagnostics output", 
            default_value = 1.0 | nbody_system.time
        )

        self.stopping_conditions.define_parameters(object)
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        
        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration, nbody_system.acceleration, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )

        object.add_method(
            "get_epsilon_squared",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )

        object.add_method(
            "set_epsilon_squared",
            (nbody_system.length * nbody_system.length, ),
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

        object.add_method(
            "get_theta_for_tree",
            (),
            (units.none, object.ERROR_CODE,)
        )

        object.add_method(
            "set_theta_for_tree",
            (units.none, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_use_self_gravity",
            (),
            (units.none, object.ERROR_CODE,)
        )

        object.add_method(
            "set_use_self_gravity",
            (units.none, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_ncrit_for_tree",
            (),
            (units.none, object.ERROR_CODE,)
        )

        object.add_method(
            "set_ncrit_for_tree",
            (units.none, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_dt_dia",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_dt_dia",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )

        self.stopping_conditions.define_methods(object)
        
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_method('EDIT', 'post_init_setters')

    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)

    
