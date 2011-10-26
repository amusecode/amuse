import os.path
from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics

class BonsaiInterface(CodeInterface, GravitationalDynamicsInterface):
    include_headers = ['worker_code.h']
    
    def name_of_worker(self,mode):
        return 'bonsai_worker'
    
    def __init__(self, mode=None, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_worker(mode), **options)
        self.set_src_directory(os.path.join(os.path.dirname(__file__), 'src', ''))
    
    @legacy_function
    def set_src_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('src_directory', dtype='string', direction=function.IN,
            description = "The path to the Bonsai src directory.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
    
    @legacy_function
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The current timestep for the system")
        function.result_type = 'int32'
        return function


class Bonsai(GravitationalDynamics):
    
    def __init__(self, unit_converter = None, **options):
        legacy_interface = BonsaiInterface(**options)
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            unit_converter,
            **options
        )
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2", 
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.05**2 | nbody_system.length**2
        )
        
        object.add_method_parameter(
            "get_time_step", 
            "set_time_step", 
            "timestep", 
            "timestep for the system", 
            default_value = 1.0 / 64 | nbody_system.time
        )
    
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        object.add_method(
            "set_time",
            (nbody_system.time,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "set_time_step",
            (nbody_system.time,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length**2, object.ERROR_CODE,)
        )
        object.add_method(
            "set_eps2",
            (nbody_system.length**2,),
            (object.ERROR_CODE,)
        )
    

