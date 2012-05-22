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
            "set_tolerance",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tolerance",
            ( ),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
