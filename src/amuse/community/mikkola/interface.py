from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from numpy import pi
from amuse.units import constants

class MikkolaInterface(CodeInterface,
                       GravitationalDynamicsInterface):
    
#    include_headers = ['worker_code.h']
    use_modules = ['Mikkola',]
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="mikkola_worker", **keyword_arguments)
    
    @legacy_function
    def set_time_step():
        """
        Set the model timestep.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The current model timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time step was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function
        
    @legacy_function
    def get_lightspeed():
        """
        Get the lightspeed value, the lightspeed scales the units (like G=1)
        and limits the valid velocity terms
        """
        function = LegacyFunctionSpecification()
        function.addParameter('lightspeed', dtype='float64',
                              direction=function.OUT,
            description = "value for the lightspeed")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_lightspeed():
        """
        Set the lightspeed value, the lightspeed scales the units (like G=1)
        and limits the valid velocity terms
        """
        function = LegacyFunctionSpecification()
        function.addParameter('lightspeed', dtype='float64',
                              direction=function.IN,
            description = "value for the lightspeed")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
        
    @legacy_function
    def get_tolerance():
        """
        Retrieve the accurancy parameter for the evolve
        """
        function = LegacyFunctionSpecification()
        function.addParameter('tolerance', dtype='float64',
                              direction=function.OUT,
            description = "tolerance for the evolve")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_tolerance():
        """
        Set the accurancy parameter for the evolve
        """
        function = LegacyFunctionSpecification()
        function.addParameter('lightspeed', dtype='float64',
                              direction=function.IN,
            description = "tolerance for the evolve")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    
    @legacy_function
    def get_radiated_gravitational_energy():
        """
        Retrieve the current radiated gravitational energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('radiated_gravitational_energy', dtype='float64', direction=function.OUT,
            description = "The energy radiated by gravitational waves")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the energy was set
        -1 - ERROR
            Energy could not be provided
        """
        return function
    
    
    @legacy_function
    def get_total_energy():
        """
        Retrieve the current total energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('energy', dtype='float64', direction=function.OUT,
            description = "The energy radiated by gravitational waves")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the energy was set
        -1 - ERROR
            Energy could not be provided
        """
        return function
        
class Mikkola(GravitationalDynamics):

    def __init__(self, convert_nbody=None, **options):
        GravitationalDynamics.__init__(
            self, 
            MikkolaInterface(**options),
            convert_nbody,
            **options
        )
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        
        if not self.unit_converter is None:
            value=self.unit_converter.to_nbody(constants.c)
            self.parameters._original.lightspeed = value

        return result
        
    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)
        
        object.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "timestep",
            "initial timestep for iteration", 
            default_value = 1.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_lightspeed", 
            "set_lightspeed",
            "lightspeed", 
            "lightspeed used in the code", 
            default_value = 1.0 | nbody_system.length / nbody_system.time
        )
        object.add_method_parameter(
            "get_tolerance", 
            "set_tolerance",
            "tolerance", 
            "tolerance used in the code", 
            default_value = 1e13
        )
        
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        
        object.add_method(
            "set_lightspeed",
            ( nbody_system.length / nbody_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_lightspeed",
            ( ),
            (nbody_system.length / nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_time_step",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "set_tolerance",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tolerance",
            ( ),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_radiated_gravitational_energy",
            (),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, object.ERROR_CODE,)
        )


        object.add_method(
            "get_total_energy",
            (),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, object.ERROR_CODE,)
        )

    def define_properties(self, object):
        
        GravitationalDynamics.define_properties(self, object)
        object.add_property("get_radiated_gravitational_energy")

