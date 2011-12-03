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
        
    @legacy_function
    def set_mass():
        """
        Update the mass of a particle. Mass is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """        

        return function     
        
    @legacy_function
    def set_state():
        """
        Update the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, position and velocity) is updated.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle")  
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        """
        return function
  
        
    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
        return function       
        
    @legacy_function
    def set_velocity():
        """
        Set the velocity vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The current x component of the velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The current y component of the velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The current z component of the velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function        
        
    @legacy_function
    def set_acceleration():
        """
        Set the velocity vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('ax', dtype='float64', direction=function.IN, description = "The current x component of the velocity vector of the particle")
        function.addParameter('ay', dtype='float64', direction=function.IN, description = "The current y component of the velocity vector of the particle")
        function.addParameter('az', dtype='float64', direction=function.IN, description = "The current z component of the velocity vector of the particle")
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
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
    

