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
        if HAS_MPMATH:
            self.adjust_prec()

    def adjust_prec(self):
        if not HAS_MPMATH:
            raise Exception("mpmath not available")
        len_ = self.parameters.word_length
        if (len_ > mpmath.mp.prec):
            mpmath.mp.prec=len_

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
