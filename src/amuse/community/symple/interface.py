from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import SinglePointGravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class sympleInterface(CodeInterface,
                   LiteratureReferencesMixIn,
                   GravitationalDynamicsInterface,
                   StoppingConditionInterface,
                   SinglePointGravityFieldInterface):
    """
    N-body integration module with shared but variable time step
    (the same for all particles but with size changing in time),
    using the symple integration scheme.  Code is symplectic for
    fixed time steps.

    .. [#] McMillan, S., 2017
    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="symple_worker",
                                 **options)
        LiteratureReferencesMixIn.__init__(self)

    def reinitialize_particles(self):
        self.recommit_particles()

    @legacy_function
    def get_dmdt():
        """
        Retrieve the mass loss rate of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32',
                              direction=function.IN,
                              description = "Index of the particle")
        function.addParameter('dmdt', dtype='float64',
                              direction=function.OUT,
                              description = "Current mdot of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK - particle was found in the model and the information was retrieved
        -1 - ERROR - particle could not be found
        """
        return function

    @legacy_function
    def set_dmdt():
        """
        Set the mass loss rate of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32',
                              direction=function.IN,
                              description = "Index of the particle")
        function.addParameter('dmdt', dtype='float64',
                              direction=function.IN,
                              description = "New mdot of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK - particle was found in the model and the information was set
        -1 - ERROR - particle could not be found
        """
        return function

    @legacy_function
    def get_timestep():
        """
        Get the fixed timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64',
                              direction=function.OUT,
            description = "fixed timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_timestep():
        """
        Set the fixed timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64',
                              direction=function.IN,
            description = "fixed timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
        
    @legacy_function
    def get_integrator():
        """
        Get the integration scheme.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('integrator', dtype='int32',
                              direction=function.OUT,
            description = "the integration scheme")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_integrator():
        """
        Set the integration scheme.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('integrator', dtype='int32',
                              direction=function.IN,
            description = "the integration scheme")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
        
    @legacy_function
    def get_timestep_parameter():
        """
        Get the time step scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('timestep_parameter', dtype='float64',
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
    def set_timestep_parameter():
        """
        Set the time step scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('timestep_parameter', dtype='float64',
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
    
class sympleDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__


class symple(GravitationalDynamics, GravityFieldCode):

    __doc__ = sympleDoc()
    __interface__ = sympleInterface

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
            "get_integrator",
            "set_integrator", 
            "integrator", 
            "integrator for gravity calculations", 
            default_value = 2
        )
        handler.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        handler.add_method_parameter(
            "get_timestep",
            "set_timestep",
            "timestep",
            "fixed timestep", 
            default_value = 0.01 | nbody_system.time
        )
        handler.add_method_parameter(
            "get_timestep_parameter",
            "set_timestep_parameter",
            "timestep_parameter",
            "timestep scaling parameter", 
            default_value = 0.05
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
            "set_integrator",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_timestep_parameter",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_time",
            (nbody_system.time,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_timestep",
            (nbody_system.time,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_dmdt",
            (handler.NO_UNIT,
             nbody_system.mass / nbody_system.time,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_integrator",
            (),
            (handler.NO_UNIT,
             handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length,
             handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_timestep_parameter",
            (),
            (handler.NO_UNIT,
             handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_time",
            (),
            (nbody_system.time,
             handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_timestep",
            (),
            (nbody_system.time,
             handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_dmdt",
            (handler.NO_UNIT,),
            (nbody_system.mass / nbody_system.time,
             handler.ERROR_CODE,)
        )

        self.stopping_conditions.define_methods(handler)
    
    def define_particle_sets(self, handler):
        GravitationalDynamics.define_particle_sets(self, handler)

        handler.add_setter('particles', 'set_dmdt')
        handler.add_getter('particles', 'get_dmdt')

        self.stopping_conditions.define_particle_set(handler)


Symple = symple
