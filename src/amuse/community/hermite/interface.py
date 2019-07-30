from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class HermiteInterface(CodeInterface,
                       LiteratureReferencesMixIn,
                       GravitationalDynamicsInterface,
                       StoppingConditionInterface,
                       SinglePointGravityFieldInterface):
    """
    N-body integration module with shared but variable time step
    (the same for all particles but its size changing in time),
    using the Hermite integration scheme.


    .. [#] Hut, P., Makino, J. & McMillan, S., *Astrophysical Journal Letters* , **443**, L93-L96 (1995)
    """
    include_headers = ['worker_code.h', 'stopcond.h']
    __so_module__ = 'hermite0_cython'

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="hermite0_worker",
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
        function.addParameter('dt_param', dtype='float64',
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
        function.addParameter('dt_param', dtype='float64',
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
    def get_is_time_reversed_allowed():
        """
        If time reversed is allowed for this code, the code will
        calculate backwards in time if the endtime given in 
        evolve_model is less than the system time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_is_time_reversed_allowed():
        """
        If time reversed is allowed for this code, the code will
        calculate backwards in time if the endtime given in 
        evolve_model is less than the system time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.IN)
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
    def get_end_time_accuracy_factor():
        """
        Get the end time accuracy factor:
            < 0 will stop between factor * dt befor the end time and the end time
              0 will stop at exactly the end time
            > 0 will stop between the end time and factor * dt after the end time
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',
                              direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_end_time_accuracy_factor():
        """
        Set the end time accuracy factor
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',
                              direction=function.IN)
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


class Hermite(GravitationalDynamics, GravityFieldCode):

    __doc__ = HermiteDoc()
    __interface__ = HermiteInterface

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = self.__interface__(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)
        self.stopping_conditions.define_state(handler)
        

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        handler.add_method_parameter(
            "get_dt_param",
            "set_dt_param",
            "dt_param",
            "timestep scaling factor", 
            default_value = 0.03
        )
        handler.add_method_parameter(
            "get_end_time_accuracy_factor",
            "set_end_time_accuracy_factor",
            "end_time_accuracy_factor",
            """
            Get the end time accuracy factor:
                < 0.0  will stop on or before the end time (larger factor, more before)
                  0.0  will stop at exactly the end time
                > 0.0  will stop on or after the end time
                
            Valid factors are between -1.0 and 1.0
            """,
            default_value = 1.0
        )
        handler.add_method_parameter(
            "get_dt_dia",
            "set_dt_dia",
            "dt_dia", 
            "time interval between diagnostics output", 
            default_value = 1.0 | nbody_system.time
        )
        handler.add_method_parameter(
            "get_is_time_reversed_allowed",
            "set_is_time_reversed_allowed",
            "is_time_reversed_allowed", 
            "if True will calculate back in time when evolve_model end time is less than systemtime", 
            default_value = False
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
            "get_dt_param",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_dt_param",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_dt_dia",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_dt_dia",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_time",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_time",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_pair_detect_factor",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_pair_detect_factor",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        self.stopping_conditions.define_methods(handler)
    
    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)
        
        self.stopping_conditions.define_particle_set(handler)
