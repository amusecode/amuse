from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface,GravityFieldInterface
from amuse.community.interface.gd import GravitationalDynamics,GravityFieldCode

class HuaynoInterface(CodeInterface,
                      LiteratureReferencesMixIn,
                      GravitationalDynamicsInterface,
                      StoppingConditionInterface,
                      GravityFieldInterface):
    """
    HUAYNO is a code to solve the astrophysical N-body problem. It uses
    recursive Hamiltonian splitting to generate multiple-timestep integrators
    which conserve momentum to machine precision. A number of different
    integrators are available. The code has been developed within the
    AMUSE environment. It can make use of GPUs - for this an OpenCL
    version can be compiled.

    .. [#] Pelupessy, Federico I.; J\"anes, J\"urgen; Portegies Zwart, Simon, New Astronomy, Volume 17, Issue 8, p. 711-719 [2012NewA...17..711P]
    .. [#] J\"anes, J\"urgen; Pelupessy, Federico I.; Portegies Zwart, Simon, A&A, Volume 570, October 2014 (for CC, OK methods) [2014A&A...570A..20J]
    """
    include_headers = ['worker_code.h']

    MODE_OPENCL='opencl'
    MODE_OPENMP='openmp'

    def name_of_worker(self,mode):
        if mode==self.MODE_OPENCL:
            return 'huayno_worker_cl'
        if mode==self.MODE_OPENMP:
            return 'huayno_worker_mp'
        return 'huayno_worker'

    def __init__(self, mode=None, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('kinetic_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_potential_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('time_param', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('time_param', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_verbosity_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('verbosity', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_verbosity_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('verbosity', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_inttype_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('inttype', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_inttype_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('inttype', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_eps2_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_eps2_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_accel_zero_mass_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('accelerate_zero_mass', dtype='b', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_accel_zero_mass_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('accelerate_zero_mass', dtype='b', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_opencl_device_type():
        function = LegacyFunctionSpecification()
        function.addParameter('opencl_device_type', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_opencl_device_type():
        function = LegacyFunctionSpecification()
        function.addParameter('opencl_device_type', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function


    def set_eps2(self, e):
        return self.set_eps2_parameter(e)

    def get_eps2(self):
        return self.get_eps2_parameter()

    @legacy_function
    def get_evolve_statistics():
        function = LegacyFunctionSpecification()
        function.addParameter('ttot', dtype='int64', direction=function.OUT)
        function.addParameter('ktot', dtype='int64', direction=function.OUT)
        function.addParameter('dtot', dtype='int64', direction=function.OUT)
        function.addParameter('tstot', dtype='int64', direction=function.OUT)
        function.addParameter('kstot', dtype='int64', direction=function.OUT)
        function.addParameter('dstot', dtype='int64', direction=function.OUT)
        function.result_type = 'i'
        return function

class Huayno(GravitationalDynamics,GravityFieldCode):

    __interface__ = HuaynoInterface

    class inttypes(object):
        """
        CONSTANT# = constant global timestep, of different order
        SHARED# = shared, but varying global timestep, of different order
        SHARED#_COLLISION = shared, but varying global timestep, of different order with collision detection
        CC_.. = various variant of connected component (termination with KEPLER or Bulirsch-stoer, see paper Janes
        CCC_... = with centering of subsys
        OK = Optimal Kick (see paper Janes)
        PASS, HOLD, BRIDGE and variants= momentum conserving individual timestepping see paper Pelupessy
        NAIVE = naive implementation of individual timestepping
        others are experimental, testing, development 
        """
        @classmethod
        def _list(cls):
              return set([x for x in cls.__dict__.keys() if not x.startswith('_')])

    all_inttypes=dict(CONSTANT = 0, SHARED2 = 1, PASS_KDK = 2, HOLD_KDK = 3, BRIDGE_KDK = 4, 
      EXTRAPOLATE = 5, PASS_DKD = 7, HOLD_DKD = 8, PPASS_DKD = 9, BRIDGE_DKD = 10,
      CC = 11, CC_KEPLER = 12, OK = 13, KEPLER = 14, SHARED4 = 15, FOURTH_M4 = 16, FOURTH_M5 = 17,
      SHARED6 = 18, SHARED8 = 19, SHARED10 = 20, SHAREDBS = 21, CCC = 22, CCC_KEPLER = 23,
      CC_BS = 24, CCC_BS = 25, BS_CC_KEPLER = 26, CC_BSA = 27, CCC_BSA = 28, SHARED2_COLLISIONS = 29,
      SHARED4_COLLISIONS = 30, SHARED6_COLLISIONS = 31, SHARED8_COLLISIONS = 32, 
      SHARED10_COLLISIONS = 33, CONSTANT2 = 34, CONSTANT4 = 35, CONSTANT6 = 36, 
      CONSTANT8 = 37, CONSTANT10 = 38, ERROR_CONTROL=39, CC_SHARED10=40, CCC_SHARED10=41)

    for key, val in all_inttypes.items():
      setattr(inttypes, key, val)

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)
        legacy_interface = self.__interface__(**options)
#        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def set_integrator(self, name):
        return self.set_inttype_parameter(self.all_inttypes[name])
    
    def get_integrator(self):
        value= self.get_inttype_parameter()
        for key, index in self.all_inttypes.items():
            if value == index:
                return key
        return "unknown"

    def define_parameters(self, handler):

        self.stopping_conditions.define_parameters(handler)

        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared",
            "smoothing parameter for gravity calculations",
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )

        handler.add_method_parameter(
            "get_timestep_parameter",
            "set_timestep_parameter",
            "timestep_parameter",
            "timestep parameter for gravity calculations",
            default_value = 0.03
        )

        handler.add_method_parameter(
            "get_timestep",
            "set_timestep",
            "timestep",
            "timestep for evolve calls",
            default_value = 0.0 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_verbosity_parameter",
            "set_verbosity_parameter",
            "verbosity_parameter",
            "verbosity parameter (0 mean silent)",
            default_value = 0
        )

        handler.add_boolean_parameter(
            "get_accel_zero_mass_parameter",
            "set_accel_zero_mass_parameter",
            "accelerate_zero_mass",
            "accelerate zero mass particle interactions (should always be true, except for testing)",
            default_value = True
        )

        inttypes=sorted([(getattr(self.inttypes,t),t ) 
                   for i,t in enumerate(sorted(self.inttypes._list()))])

        handler.add_method_parameter(
            "get_inttype_parameter",
            "set_inttype_parameter",
            "inttype_parameter",
            "integrator method to use, this can be one of: "+
             ",".join( ["{0}={1}".format(i, t) for i,t in inttypes]),
            #~ default_value = 8
        )

        handler.add_method_parameter(
            "get_integrator",
            "set_integrator",
            "integrator",
            "integrator method to use, this can be one of: "+
             ",".join( ["{0}".format(t) for i,t in inttypes]),
             #~ default_value="HOLD_DKD"
        )

        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

        handler.add_method_parameter(
            "get_opencl_device_type",
            "set_opencl_device_type",
            "opencl_device_type",
            "set preferred OpenCL device type (0=default, 1=cpu, 2=gpu)",
            default_value = 0
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_timestep_parameter",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_timestep_parameter",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_timestep",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_timestep",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_inttype_parameter",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_inttype_parameter",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        self.stopping_conditions.define_methods(handler)

    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        self.stopping_conditions.define_particle_set(handler)

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)

        handler.add_method('RUN', 'get_kinetic_energy')
        handler.add_method('RUN', 'get_potential_energy')

        self.stopping_conditions.define_state(handler)
