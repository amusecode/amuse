from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

class AdaptbInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, 
        StoppingConditionInterface, CodeWithDataDirectories):
    """
    Adaptb (Accurate Dynamics with Arbitrary Precision by Tjarda Boekholt)
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="adaptb_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
    
    @legacy_function
    def get_adaptb_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('adaptb_output_directory', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_adaptb_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('adaptb_output_directory', dtype='string', direction=function.IN)
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

    @legacy_function
    def new_particle_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('identity_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='string', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='string', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='string', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='string', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='string', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='string', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='string', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('radius', dtype='string', direction=function.IN, description = "The radius of the particle", default='0')
        function.result_type = 'int32'
        return function

    def new_particle(self, mass, x,y,z, vx,vy,vz, radius = 0):
        if isinstance(mass, str):
            return self.new_particle_string(mass, x,y,z, vx,vy,vz, radius = str(radius))
        else:
            return self.new_particle_float64(mass, x,y,z, vx,vy,vz, radius = radius)

    @legacy_function
    def get_bs_tolerance_string():
        function = LegacyFunctionSpecification()
        function.addParameter('bs_tolerance', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_bs_tolerance_string():
        function = LegacyFunctionSpecification()
        function.addParameter('bs_tolerance', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_bs_tolerance_float64():
        function = LegacyFunctionSpecification()
        function.addParameter('bs_tolerance', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_bs_tolerance_float64():
        function = LegacyFunctionSpecification()
        function.addParameter('bs_tolerance', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_word_length():
        function = LegacyFunctionSpecification()
        function.addParameter('word_length', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_word_length():
        function = LegacyFunctionSpecification()
        function.addParameter('word_length', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dt_print():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_print', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_dt_print():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_print', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_max_cpu_time():
        function = LegacyFunctionSpecification()
        function.addParameter('max_cpu_time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_max_cpu_time():
        function = LegacyFunctionSpecification()
        function.addParameter('max_cpu_time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    


class Adaptb(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = AdaptbInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        ensure_data_directory_exists(self.output_directory)
        self.parameters.adaptb_output_directory = self.output_directory
        return result
    
    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        self.stopping_conditions.define_parameters(object)
        
        object.add_method_parameter(
            "get_bs_tolerance_float64", 
            "set_bs_tolerance_float64",
            "bs_tolerance", 
            "Error tolerance of the Bulirsch-Stoer integrator", 
            default_value = 1.0e-6
        )
        
        object.add_method_parameter(
            "get_eps2", 
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations, usage is not recommended for Adaptb", 
            default_value = 0.0 | nbody_system.length**2
        )
        
        object.add_method_parameter(
            "get_dt_print", 
            "set_dt_print",
            "dt_print", 
            "dt_print, regular print interval to show status (% complete) of evolve_model", 
            default_value = 0.1 | nbody_system.time
        )
    
        object.add_method_parameter(
            "get_word_length", 
            "set_word_length",
            "word_length", 
            "The word length, or number of bits, used for the arbitrary precision calculations", 
            default_value = 64
        )
        
        object.add_method_parameter(
            "get_adaptb_output_directory", 
            "set_adaptb_output_directory",
            "adaptb_output_directory", 
            "Path to the directory where Adaptb stores its output", 
            default_value = "./"
        )
        
        object.add_method_parameter(
            "get_max_cpu_time", 
            "set_max_cpu_time",
            "time_limit_cpu", 
            "The cpu-time limit, the maximum amount of time Adaptb is allowed to run for.", 
            default_value = 3600.0 | units.s
        )
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
        
        object.add_method("get_bs_tolerance_float64", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_bs_tolerance_float64", (object.NO_UNIT, ), (object.ERROR_CODE,))
        
        object.add_method("get_eps2", (), (nbody_system.length**2, object.ERROR_CODE,))
        object.add_method("set_eps2", (nbody_system.length**2, ), (object.ERROR_CODE,))
    
        object.add_method("get_dt_print", (), (nbody_system.time, object.ERROR_CODE,))
        object.add_method("set_dt_print", (nbody_system.time, ), (object.ERROR_CODE,))
    
        object.add_method("get_word_length", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_word_length", (object.NO_UNIT, ), (object.ERROR_CODE,))
        
        object.add_method("get_adaptb_output_directory", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_adaptb_output_directory", (object.NO_UNIT, ), (object.ERROR_CODE,))
        
        object.add_method("get_max_cpu_time", (), (units.s, object.ERROR_CODE,))
        object.add_method("set_max_cpu_time", (units.s, ), (object.ERROR_CODE,))
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)


