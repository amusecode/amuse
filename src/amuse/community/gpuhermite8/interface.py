from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
from amuse.community.interface.gd import GravityFieldCode
try:
    import mpmath
    HAS_MPMATH=True
except ImportError:
    HAS_MPMATH=False

class gpuhermite8Interface(CodeInterface, GravitationalDynamicsInterface, StoppingConditionInterface):

    include_headers = ['worker_code.h', 'stopcond.h']
    def __init__(self, **keyword_arguments):
#        print("gpuhermite8Interface")
        CodeInterface.__init__(self, name_of_the_worker="gpuhermite8_worker", **keyword_arguments)

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
    @legacy_function
    def get_const_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('const_time_step', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_const_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('const_time_step', dtype='float64', direction=function.IN)
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
    def get_potential_energy_string():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy_string', dtype='string', direction=function.OUT)
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
    def get_total_mass_string():
        function = LegacyFunctionSpecification()
        function.addParameter('M', dtype='string', direction=function.OUT)
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

    @legacy_function
    def add_step_acceleration_float64():
        """
         Add acceleration to a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('a_step_x', dtype='float64', direction=function.IN, description = "Add acceleration vector to the particle (valid for next evaluation step only)")
        function.addParameter('a_step_y', dtype='float64', direction=function.IN, description = "Add acceleration vector to the particle (valid for next evaluation step only)")
        function.addParameter('a_step_z', dtype='float64', direction=function.IN, description = "Add acceleration vector to the particle (valid for next evaluation step only)")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function


    def new_particle(self, mass, x,y,z, vx,vy,vz, radius = 0):
        if isinstance(mass, str):
            return self.new_particle_string(mass, x,y,z, vx,vy,vz, radius = str(radius))
        else:
            return self.new_particle_float64(mass, x,y,z, vx,vy,vz, radius = radius)



class Gpuhermite8(GravitationalDynamics, GravityFieldCode):
    def __init__(self, convert_nbody = None, **options):
#        print("Gpuhermite8")
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = gpuhermite8Interface(**options)
        self.legacy_doc = legacy_interface.__doc__
#        InCodeComponentImplementation.__init__(self,  gpuhermite8Interface(**options), **options)

class gpuhermite8(GravitationalDynamics, GravityFieldCode):
    def __init__(self, convert_nbody = None, **options):
#        print("gpuhermite8")
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = gpuhermite8Interface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )
        self.convert_nbody=convert_nbody
        self.not_adjusted = True
        if HAS_MPMATH:
            self.adjust_prec()

    def adjust_prec(self):
        if self.not_adjusted:
            if not HAS_MPMATH:
                raise Exception("mpmath not available")
            len_ = self.parameters.word_length
            if (len_ > mpmath.mp.prec):
                mpmath.mp.prec=len_
            self.not_adjusted = False

    def get_potential_energy_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.energy).number
#        print (a,self.convert_nbody.to_si(nbody_system.energy))
        b=mpmath.mpf(self.get_potential_energy_string())*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit

    def get_total_energy_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.energy).number
#        print (a)
        c= self.get_total_energy_string()
        b=mpmath.mpf(c)*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit

    def get_kinetic_energy_p_si(self):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.energy).number
#        print (a)
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
        b=mpmath.mpf(self.get_mass_string(index+1))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit

    def get_radius_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.length).number
        b=mpmath.mpf(self.get_radius_string(index+1))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit

    def get_velocity_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.speed).number
        b=mpmath.matrix(self.get_velocity_string(index+1))*a
        return b #* self.convert_nbody.to_si(nbody_system.energy).unit

    def get_position_p_si(self,index):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.length).number
        b=mpmath.matrix(self.get_position_string(index+1))*a
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

    def add_step_acceleration(self, index, a_step_x, a_step_y, a_step_z):
        self.adjust_prec()
        a=self.convert_nbody.to_si(nbody_system.acceleration)#.number
        return self.add_step_acceleration_float64(index+1 , a_step_x/a, a_step_y/a, a_step_z/a)

    def set_state_p_si(self,index,m,x,y,z,vx,vy,vz,radius):
        self.adjust_prec()
        um=self.convert_nbody.to_si(nbody_system.mass).number
        ul=self.convert_nbody.to_si(nbody_system.length).number
        us=self.convert_nbody.to_si(nbody_system.speed).number
        return self.set_state_string(index+1 ,str(m/um) ,str(x/ul) ,
                                        str(y/ul) ,str(z/ul) ,
                                        str(vx/us) ,str(vy/us) ,
                                        str(vz/us) ,str(radius/ul))

    def set_mass_p_si(self,index,m):
        self.adjust_prec()
        um=self.convert_nbody.to_si(nbody_system.mass).number
        self.set_mass_string(index+1 ,str(m/um))

    def set_radius_p_si(self,index,radius):
        self.adjust_prec()
        ul=self.convert_nbody.to_si(nbody_system.length).number
        self.set_radius_string(index+1 ,str(radius/ul))


    def set_position_p_si(self,index,x,y,z):
        self.adjust_prec()
        ul=self.convert_nbody.to_si(nbody_system.length).number
        self.set_position_string(index+1 ,str(x/ul) ,str(y/ul) ,str(z/ul))

    def set_velocity_p_si(self,index,vx,vy,vz):
        self.adjust_prec()
        us=self.convert_nbody.to_si(nbody_system.speed).number
        self.set_velocity_string(index+1 ,str(vx/us) ,str(vy/us) ,str(vz/us))


    def get_state_p_si(self,index):
        self.adjust_prec()
        b=mpmath.matrix(self.get_state_string(index+1))
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
            "get_word_length",
            "set_word_length",
            "word_length",
            "The word length, or number of bits for the mantissa, used for the arbitrary precision calculations (#digits = log10(2**# bits) ",
            default_value = 53
        )

        handler.add_method_parameter(
            "get_const_time_step",
            "set_const_time_step",
            "const_timestep",
            "constant timestep for iteration",
            default_value = 0.7 | nbody_system.time
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        self.stopping_conditions.define_methods(handler)

        handler.add_method("get_word_length", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_word_length", (handler.NO_UNIT, ), (handler.ERROR_CODE,))

        handler.add_method("get_const_time_step", (), (nbody_system.time , handler.ERROR_CODE,))
        handler.add_method("set_const_time_step", ( nbody_system.time ), (handler.ERROR_CODE,))

    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)

        handler.add_method('INITIALIZED', 'get_word_length')
        handler.add_method('EDIT', 'get_word_length')
        handler.add_method('UPDATE', 'get_word_length')
        handler.add_method('RUN', 'get_word_length')

        handler.add_method('RUN', 'get_total_energy_string')
        handler.add_method('RUN', 'get_kinetic_energy_string')
        handler.add_method('RUN', 'get_potential_energy_string')

        handler.add_method('RUN', 'get_state_string')
        handler.add_method('RUN', 'get_position_string')
        handler.add_method('RUN', 'get_velocity_string')
        handler.add_method('RUN', 'get_radius_string')
        handler.add_method('RUN', 'get_mass_string')
        handler.add_method('RUN', 'get_total_mass_string')
        handler.add_method('EDIT', 'new_particle_string')
        handler.add_method('UPDATE', 'new_particle_string')
        handler.add_method('EDIT', 'add_step_acceleration_float64')
        handler.add_method('UPDATE', 'add_step_acceleration_float64')

        handler.add_transition('RUN', 'UPDATE', 'new_particle_string', False)
        handler.add_transition('RUN', 'UPDATE', 'set_velocity_string', False)
        handler.add_transition('RUN', 'UPDATE', 'set_position_string', False)
        handler.add_transition('RUN', 'UPDATE', 'set_radius_string', False)
        handler.add_transition('RUN', 'UPDATE', 'set_mass_string', False)
        handler.add_transition('RUN', 'UPDATE', 'set_state_string', False)
        handler.add_transition('RUN', 'UPDATE', 'add_step_acceleration_float64', False)

