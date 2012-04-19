from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics

class HermiteInterface(CodeInterface,
                       LiteratureReferencesMixIn,
                       GravitationalDynamicsInterface,
                       StoppingConditionInterface):
    """
    N-body integration module with shared but variable time step
    (the same for all particles but its size changing in time),
    using the Hermite integration scheme.


    .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """
    include_headers = ['worker_code.h', 'stopcond.h']


    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="hermite_worker",
                                 **options)
        LiteratureReferencesMixIn.__init__(self)

    def reinitialize_particles(self):
        self.recommit_particles()
    
    @legacy_function
    def get_dt_dia():
        """
        Get the time interval between diagnostics output.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.OUT,
            description = "time interval between diagnostic outputs")
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
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.IN,
            description = "the time interval between diagnostics output")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_dt_param():
        """
        Get the timestep scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.OUT,
            description = "the timestep scaling factor")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_dt_param():
        """
        Set the timestep scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.IN,
            description = "the timestep scaling factor")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_time():
        """
        Get the current simulation time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64',
                              direction=function.OUT,
            description = "the current simulation time")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_time():
        """
        Set the current simulation time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64',
                              direction=function.IN,
            description = "the current simulation time")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_pair_detect_factor():
        """
        Get pair detection sphere radius factor (units particle radius).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('pair_detect_factor', dtype='float64',
                              direction=function.OUT,
            description = "pair detection radius ")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def set_pair_detect_factor():
        """
        Set pair detection sphere radius factor (units particle radius).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('pair_detect_factor', dtype='float64',
                              direction=function.IN,
            description = "pair detection radius ")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    

class HermiteDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__


class Hermite(GravitationalDynamics):

    __doc__ = HermiteDoc()

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = HermiteInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_dt_param",
            "set_dt_param",
            "dt_param",
            "timestep scaling factor", 
            default_value = 0.03
        )
        object.add_method_parameter(
            "get_dt_dia",
            "set_dt_dia",
            "dt_dia", 
            "time interval between diagnostics output", 
            default_value = 1.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_time",
            "set_time",
            "time",
            "current simulation time", 
            default_value = 0.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_pair_detect_factor",
            "set_pair_detect_factor",
            "pair_factor",
            "radius factor for pair detection", 
            default_value = 1.0
        )

        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length,
             nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )

        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length,
             nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration,
             nbody_system.acceleration, object.ERROR_CODE)
        )
        
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_dt_param",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_dt_param",
            (object.NO_UNIT, ),
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
        
        object.add_method(
            "get_time",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_pair_detect_factor",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_pair_detect_factor",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        self.stopping_conditions.define_methods(object)
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        
        self.stopping_conditions.define_particle_set(object)
