from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

class AdaptbInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, StoppingConditionInterface):
    """
    Adaptb (Accurate Dynamics with Arbitrary Precision by Tjarda Boekholt)
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="adaptb_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
    
    @option(type="string", sections=('data',))
    def output_data_root_directory(self):
        """
        The root directory of the output data,
        read - write directory
        """
        return os.path.join(get_amuse_root_dir(), 'data')
        
    @option(type="string")
    def output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'adaptb', 'output')

    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return self.output_directory
  
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
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('identity_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('mass', dtype='string', direction=function.IN, description = "The mass of the particle")
        function.addParameter('radius', dtype='string', direction=function.IN, description = "The radius of the particle")
        function.addParameter('x', dtype='string', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='string', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='string', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='string', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='string', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='string', direction=function.IN, description = "The initial velocity vector of the particle")
        function.result_type = 'int32'
        return function

    def new_particle(self, mass, radius, x,y,z, vx,vy,vz):
        if isinstance(mass, str):
            return self.new_particle_string(mass, radius, x,y,z, vx,vy,vz)
        else:
            return self.new_particle_float64(mass, radius, x,y,z, vx,vy,vz)

    @legacy_function
    def get_bs_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('bs_tolerance', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_bs_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('bs_tolerance', dtype='string', direction=function.IN)
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
    def get_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_dt_print():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_print', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_dt_print():
        function = LegacyFunctionSpecification()
        function.addParameter('dt_print', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_max_cpu_time():
        function = LegacyFunctionSpecification()
        function.addParameter('max_cpu_time', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_max_cpu_time():
        function = LegacyFunctionSpecification()
        function.addParameter('max_cpu_time', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    #@legacy_function
    #def evolve_model():
        #function = LegacyFunctionSpecification()
        #function.addParameter('t', dtype='string', direction=function.IN)
        #function.result_type = 'int32'
        #return function

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

    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        self.stopping_conditions.define_parameters(object)
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        self.stopping_conditions.define_methods(object)
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object, 'particles')


