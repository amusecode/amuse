from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

"""
currently setting the particle (and possibly model time) as strings (ie to conserve 
precision) is not yet supported fully (no high level, low level untested)
"""

class BrutusInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, 
        StoppingConditionInterface, CodeWithDataDirectories):
    """
    Brutus (Brute force N-body code)
        .. [#] Boekholt, Tjarda and Portegies Zwart, Simon,On the reliability of N-body simulations, Computational Astrophysics and Cosmology, Volume 2, article id.2, 21 pp.
    
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    ####
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="brutus_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    ####
    @legacy_function
    def get_brutus_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('brutus_output_directory', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_brutus_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('brutus_output_directory', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    ####
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

    ####
    @legacy_function
    def get_bs_tolerance_string():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_bs_tolerance_string():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_bs_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_bs_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    ####
    @legacy_function
    def get_word_length():
        function = LegacyFunctionSpecification()
        function.addParameter('numBits', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_word_length():
        function = LegacyFunctionSpecification()
        function.addParameter('numBits', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    ####
    @legacy_function
    def get_eta_string():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_param', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta_string():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_param', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_param', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_param', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    ####
    @legacy_function
    def get_t_string():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_t_string():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_t():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_t():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

class Brutus(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = BrutusInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        self.parameters.brutus_output_directory = self.output_directory
        return result
    
    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        self.stopping_conditions.define_parameters(object)
        
        object.add_method_parameter(
            "get_bs_tolerance", 
            "set_bs_tolerance",
            "bs_tolerance", 
            "Error tolerance of the Bulirsch-Stoer integrator", 
            default_value = 1.0e-8
        )

        object.add_method_parameter(
            "get_word_length", 
            "set_word_length",
            "word_length", 
            "The word length, or number of bits for the mantissa, used for the arbitrary precision calculations (#digits = log10(2**# bits) ", 
            default_value = 72
        )
                
        object.add_method_parameter(
            "get_eta", 
            "set_eta",
            "dt_param", 
            "dt_param, the time-step parameter for the adaptive time-step criterion", 
            default_value = 0.24
        )
            
        object.add_method_parameter(
            "get_brutus_output_directory", 
            "set_brutus_output_directory",
            "brutus_output_directory", 
            "Path to the directory where Brutus stores its output", 
            default_value = "./"
        )
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
        
        object.add_method("get_bs_tolerance", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_bs_tolerance", (object.NO_UNIT, ), (object.ERROR_CODE,))
 
        object.add_method("get_word_length", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_word_length", (object.NO_UNIT, ), (object.ERROR_CODE,))

        object.add_method("get_eta", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_eta", (object.NO_UNIT, ), (object.ERROR_CODE,))
        
        object.add_method("get_brutus_output_directory", (), (object.NO_UNIT, object.ERROR_CODE,))
        object.add_method("set_brutus_output_directory", (object.NO_UNIT, ), (object.ERROR_CODE,))
            
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object)


