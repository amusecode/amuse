from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class BHTreeInterface(
    CodeInterface,
    LiteratureReferencesMixIn,
    GravitationalDynamicsInterface,
    StoppingConditionInterface,
    SinglePointGravityFieldInterface):
    """
        .. [#] Barnes, J., Hut, P., A Hierarchical O(N log N) force-calculation algorithm, *Nature*, **4**, 324 (1986)   
    """
    include_headers = ['interface.h', 'worker_code.h', 'stopcond.h']
    __so_module__ = 'bhtree_cython'
    
    def __init__(self, convert_nbody = None, mode = 'cpu', **kwargs):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **kwargs)
        """
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
            
        self.convert_nbody = convert_nbody
        """
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self, mode):
        if mode == "g6":
            return 'bhtree_worker_g6'
        elif mode == "gpu":
            return 'bhtree_worker_gpu'
        else:
            return 'bhtree_worker'

       
    
    
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
            description = "the time interval between diagnostics output",
            unit = nbody_system.time)
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
            description = "the time interval between diagnostics output",
            unit = nbody_system.time)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
class BHTree(GravitationalDynamics, GravityFieldCode):

    __interface__ = BHTreeInterface
    
    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)
        
        legacy_interface = self.__interface__(**options)
        
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )     
            
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_epsilon_squared",
            "set_epsilon_squared", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.125 | nbody_system.length * nbody_system.length
        )
        handler.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "constant timestep for iteration", 
            default_value = 0.015625 | nbody_system.time
        )
        handler.add_method_parameter(
            "get_theta_for_tree",
            "set_theta_for_tree",
            "opening_angle", 
            "opening angle, theta, for building the tree: between 0 and 1", 
            default_value = 0.75
        )
        handler.add_method_parameter(
            "get_use_self_gravity",
            "set_use_self_gravity",
            "use_self_gravity", 
            "flag for usage of self gravity, 1 or 0 (true or false)", 
            default_value = 1
        )
        handler.add_method_parameter(
            "get_ncrit_for_tree",
            "set_ncrit_for_tree",
            "ncrit_for_tree", 
            "Ncrit, the maximum number of particles sharing an interaction list", 
            default_value = 12
        )
        handler.add_method_parameter(
            "get_dt_dia",
            "set_dt_dia",
            "dt_dia", 
            "time interval between diagnostics output", 
            default_value = 1.0 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )
        
        self.stopping_conditions.define_parameters(handler)
        
    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        
        handler.add_method(
            "get_time_step",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_time_step",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )



        self.stopping_conditions.define_methods(handler)
        
    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)
        
        self.stopping_conditions.define_state(handler)


    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)

    
