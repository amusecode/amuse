from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
try:
    import mpmath
    HAS_MPMATH=True
except ImportError:
    HAS_MPMATH=False


"""
currently setting the particle (and possibly model time) as strings (ie to conserve 
precision) is not yet supported fully (no high level, low level untested)
"""

class BrutusInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, 
        StoppingConditionInterface, CodeWithDataDirectories):
    """
    Brutus (Brute force N-body code)
        .. [#] Boekholt, Tjarda and Portegies Zwart, Simon, On the reliability of N-body simulations, Computational Astrophysics and Cosmology, Volume 2, article id.2, 21 pp. [2015ComAC...2....2B]
    
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
    
    @legacy_function
    def set_state_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='string', direction=function.IN, description = "new mass of the particle")
        function.addParameter('x', dtype='string', direction=function.IN, description = "new initial position vector of the particle")
        function.addParameter('y', dtype='string', direction=function.IN, description = "new initial position vector of the particle")
        function.addParameter('z', dtype='string', direction=function.IN, description = "new initial position vector of the particle")
        function.addParameter('vx', dtype='string', direction=function.IN, description = "new initial velocity vector of the particle")
        function.addParameter('vy', dtype='string', direction=function.IN, description = "new initial velocity vector of the particle")
        function.addParameter('vz', dtype='string', direction=function.IN, description = "new initial velocity vector of the particle")
        function.addParameter('radius', dtype='string', direction=function.IN, description = "new radius of the particle", default='0')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mass_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='string', direction=function.IN, description = "new mass of the particle")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_radius_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('radius', dtype='string', direction=function.IN, description = "new radius of the particle", default='0')
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_position_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='string', direction=function.IN, description = "new position vector of the particle")
        function.addParameter('y', dtype='string', direction=function.IN, description = "new position vector of the particle")
        function.addParameter('z', dtype='string', direction=function.IN, description = "new position vector of the particle")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_velocity_string():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('vx', dtype='string', direction=function.IN, description = "new velocity vector of the particle")
        function.addParameter('vy', dtype='string', direction=function.IN, description = "new velocity vector of the particle")
        function.addParameter('vz', dtype='string', direction=function.IN, description = "new velocity vector of the particle")
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
        
    @legacy_function
    def get_potential_energy_string():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy_string', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function 
        
    @legacy_function
    def get_total_mass_string():
        function = LegacyFunctionSpecification()
        function.addParameter('M', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function 
       
    @legacy_function
    def get_total_energy_string():
        function = LegacyFunctionSpecification()
        function.addParameter('total_energy_string', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function 
       
    @legacy_function
    def get_kinetic_energy_string():
        function = LegacyFunctionSpecification()
        function.addParameter('kinetic_energy_string', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_velocity_string():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('vx', dtype='string', direction=function.OUT)
        function.addParameter('vy', dtype='string', direction=function.OUT)
        function.addParameter('vz', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_position_string():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='string', direction=function.OUT)
        function.addParameter('y', dtype='string', direction=function.OUT)
        function.addParameter('z', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function  
    
    @legacy_function
    def get_mass_string():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('m', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_radius_string():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('radius', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_state_string():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('m', dtype='string', direction=function.OUT)
        function.addParameter('x', dtype='string', direction=function.OUT)
        function.addParameter('y', dtype='string', direction=function.OUT)
        function.addParameter('z', dtype='string', direction=function.OUT)
        function.addParameter('vx', dtype='string', direction=function.OUT)
        function.addParameter('vy', dtype='string', direction=function.OUT)
        function.addParameter('vz', dtype='string', direction=function.OUT)
        function.addParameter('radius', dtype='string', direction=function.OUT)
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
        self.convert_nbody=convert_nbody
        if HAS_MPMATH:
            self.adjust_prec()
            
    def adjust_prec(self):
        if not HAS_MPMATH:
            raise Exception("mpmath not available")
        len_ = self.parameters.word_length
        if (len_ > mpmath.mp.prec):
            mpmath.mp.prec=len_
        
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        self.parameters.brutus_output_directory = self.output_directory
        return result
    
    def get_potential_energy_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.energy).number
        b=mpmath.mpf(self.get_potential_energy_string())*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit
    
    def get_total_energy_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.energy).number
        b=mpmath.mpf(self.get_total_energy_string())*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit
    
    def get_kinetic_energy_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.energy).number
        b=mpmath.mpf(self.get_kinetic_energy_string())*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit
    
    def get_total_mass_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.mass).number
        b=mpmath.mpf(self.get_total_mass_string())*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit
        
    def get_mass_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.mass).number
        b=mpmath.mpf(self.get_mass_string(index))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit
        
    def get_radius_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.length).number
        b=mpmath.mpf(self.get_radius_string(index))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit
    
    def get_time_p_si(self):
        self.adjust_prec()
        return mpmath.mpf(self.get_time().number)

    def get_velocity_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.speed).number
        b=mpmath.matrix(self.get_velocity_string(index))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit   

    def get_position_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.length).number
        b=mpmath.matrix(self.get_position_string(index))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit   

    def new_particle_p_si(self,m,x,y,z,vx,vy,vz,radius):
        self.adjust_prec()
        um=self.convert_nbody.to_si(nbody_system.mass).number
        ul=self.convert_nbody.to_si(nbody_system.length).number
        us=self.convert_nbody.to_si(nbody_system.speed).number
        return self.new_particle_string(str(m/um) ,str(x/ul) ,
                                        str(y/ul) ,str(z/ul) ,
                                        str(vx/us) ,str(vy/us) ,
                                        str(vz/us) ,str(radius/ul))
  
    def set_state_p_si(self,index,m,x,y,z,vx,vy,vz,radius):
        self.adjust_prec()
        um=self.convert_nbody.to_si(nbody_system.mass).number
        ul=self.convert_nbody.to_si(nbody_system.length).number
        us=self.convert_nbody.to_si(nbody_system.speed).number
        return self.set_state_string(index ,str(m/um) ,str(x/ul) ,
                                        str(y/ul) ,str(z/ul) ,
                                        str(vx/us) ,str(vy/us) ,
                                        str(vz/us) ,str(radius/ul))
        
    def set_mass_p_si(self,index,m):
        self.adjust_prec()
        um=self.convert_nbody.to_si(nbody_system.mass).number
        self.set_mass_string(index ,str(m/um))

    def set_radius_p_si(self,index,radius):
        self.adjust_prec()
        ul=self.convert_nbody.to_si(nbody_system.length).number
        self.set_radius_string(index ,str(radius/ul))
        
        
    def set_position_p_si(self,index,x,y,z):
        self.adjust_prec()
        ul=self.convert_nbody.to_si(nbody_system.length).number
        self.set_position_string(index ,str(x/ul) ,str(y/ul) ,str(z/ul))      
        
    def set_velocity_p_si(self,index,vx,vy,vz):
        self.adjust_prec()
        us=self.convert_nbody.to_si(nbody_system.speed).number
        self.set_velocity_string(index ,str(vx/us) ,str(vy/us) ,str(vz/us))
        
        
    def get_state_p_si(self,index):
        self.adjust_prec()
        b=mpmath.matrix(self.get_state_string(index))
        a=self.convert_nbody.to_si(nbody_system.mass).number
        b[0]=b[0]*a
        a=self.convert_nbody.to_si(nbody_system.length).number
        b[1]=b[1]*a
        b[2]=b[2]*a
        b[3]=b[3]*a
        a=self.convert_nbody.to_si(nbody_system.speed).number
        b[4]=b[4]*a
        b[5]=b[5]*a
        b[6]=b[6]*a
        a=self.convert_nbody.to_si(nbody_system.length).number
        b[7]=b[7]*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit   

    def define_parameters(self, handler):
        GravitationalDynamics.define_parameters(self, handler)
        self.stopping_conditions.define_parameters(handler)
        
        handler.add_method_parameter(
            "get_bs_tolerance", 
            "set_bs_tolerance",
            "bs_tolerance", 
            "Error tolerance of the Bulirsch-Stoer integrator", 
            default_value = 1.0e-8
        )

        handler.add_method_parameter(
            "get_word_length", 
            "set_word_length",
            "word_length", 
            "The word length, or number of bits for the mantissa, used for the arbitrary precision calculations (#digits = log10(2**# bits) ", 
            default_value = 72
        )
                
        handler.add_method_parameter(
            "get_eta", 
            "set_eta",
            "dt_param", 
            "dt_param, the time-step parameter for the adaptive time-step criterion", 
            default_value = 0.24
        )
            
        handler.add_method_parameter(
            "get_brutus_output_directory", 
            "set_brutus_output_directory",
            "brutus_output_directory", 
            "Path to the directory where Brutus stores its output", 
            default_value = "./"
        )
        
    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        self.stopping_conditions.define_methods(handler)
        
        handler.add_method("get_bs_tolerance", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_bs_tolerance", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
 
        handler.add_method("get_word_length", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_word_length", (handler.NO_UNIT, ), (handler.ERROR_CODE,))

        handler.add_method("get_eta", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_eta", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_brutus_output_directory", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_brutus_output_directory", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
            
    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)


