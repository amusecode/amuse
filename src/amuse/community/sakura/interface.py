from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

MODULES_MISSING = False

class SakuraInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, 
        StoppingConditionInterface, CodeWithDataDirectories):
    """
    Sakura (Ferrari et al. 2014)
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="sakura_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    @legacy_function
    def get_sakura_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('sakura_output_directory', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_sakura_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('sakura_output_directory', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle_float64():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('identity_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle", default = 0)
        function.result_type = 'int32'
        return function

    def new_particle(self, mass, x,y,z, vx,vy,vz, radius = 0):
        return self.new_particle_float64(mass, x,y,z, vx,vy,vz, radius = radius)

    ##############################################

    @legacy_function
    def get_t_begin():
        function = LegacyFunctionSpecification()
        function.addParameter('t_begin', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_t_begin():
        function = LegacyFunctionSpecification()
        function.addParameter('t_begin', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dt():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_dt():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_t():
        function = LegacyFunctionSpecification()
        function.addParameter('t', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_t():
        function = LegacyFunctionSpecification()
        function.addParameter('t', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    ####################################################

class Sakura(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = SakuraInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        self.parameters.sakura_output_directory = self.output_directory
        return result
    
    def define_parameters(self, handler):
        GravitationalDynamics.define_parameters(self, handler)
        self.stopping_conditions.define_parameters(handler)
        
        ####################################################

        handler.add_method_parameter(
            "get_t_begin", 
            "set_t_begin",
            "begin_time", 
            "Time at start of simulation", 
            default_value = 0.0 | nbody_system.time
        )
        
        handler.add_method_parameter(
            "get_dt", 
            "set_dt",
            "timestep", 
            "Constant time-step size", 
            default_value = 1e-3 | nbody_system.time
        )
        
        handler.add_method_parameter(
            "get_t", 
            "set_t",
            "current_time", 
            "Current time", 
            default_value = 0.0 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_sakura_output_directory", 
            "set_sakura_output_directory",
            "sakura_output_directory", 
            "Output directory", 
            default_value = "./"
        )

        ####################################################

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        self.stopping_conditions.define_methods(handler)
        
        ####################################################

        handler.add_method("get_t_begin", (), (nbody_system.time, handler.ERROR_CODE,))
        handler.add_method("set_t_begin", (nbody_system.time, ), (handler.ERROR_CODE,))
    
        handler.add_method("get_dt", (), (nbody_system.time, handler.ERROR_CODE,))
        handler.add_method("set_dt", (nbody_system.time, ), (handler.ERROR_CODE,))

        handler.add_method("get_t", (), (nbody_system.time, handler.ERROR_CODE,))
        handler.add_method("set_t", (nbody_system.time, ), (handler.ERROR_CODE,))

        ####################################################
    
    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)


    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        
        self.stopping_conditions.define_state(handler)

