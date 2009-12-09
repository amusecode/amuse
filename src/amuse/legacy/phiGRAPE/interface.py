import numpy

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.interface.gd import NBodyGravitationalDynamicsBinding
from amuse.legacy.support.lit import LiteratureRefs

class PhiGRAPEInterface(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """
        .. [#] Refs not included yet
    """    
    
    MODE_G6LIB = 'g6lib'
    MODE_GPU   = 'gpu'
    MODE_GRAPE = 'grape'
    
    def __init__(self, convert_nbody = None, mode = MODE_G6LIB):
        LegacyInterface.__init__(self, name_of_the_worker = self.name_of_the_muse_worker(mode))
        LiteratureRefs.__init__(self)

    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'worker_code'
        if mode == self.MODE_GPU:
            return 'muse_worker_gpu'
        if mode == self.MODE_GRAPE:
            return 'muse_worker_grape'
            
    @legacy_function   
    def setup_module():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
        
    @legacy_function      
    def cleanup_module():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def initialize_particles():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;

    @legacy_function  
    def reinitialize_particles():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    """                
    OBSOLETE
    @legacy_function    
    def add_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    cello in gd 
    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    """    
    """
    @legacy_function   
    def get_number():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function;
    """    
    """cello in gd
    @legacy_function   
    def get_eps2():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function;
           

    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        return function

    cello in gd
    @legacy_function      
    def get_potential():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
    
        
    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        return function;

    @legacy_function      
    def get_time():
        function = LegacyFunctionSpecification()
        function.result_type = 'd'
        return function

    """
    
    @legacy_function
    def synchronize_model():
        """
        evolve all particles up to current sys time
        """
        function = LegacyFunctionSpecification() 
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            
        """
        return function  

    @legacy_function      
    def get_time_step():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function
    """
    @legacy_function      
    def set_eps():
        function = LegacyFunctionSpecification()  
        function.addParameter('eps2', dtype='d', direction=function.IN)
        return function
    """
    @legacy_function      
    def set_eta():
        function = LegacyFunctionSpecification()  
        function.addParameter('etas', dtype='d', direction=function.IN)
        function.addParameter('eta', dtype='d', direction=function.IN)
        return function

    @legacy_function      
    def set_eta_s():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='d', direction=function.IN)
        return function

    @legacy_function      
    def set_eta1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='d', direction=function.IN)
        return function


    @legacy_function      
    def get_eta():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function
        
    @legacy_function      
    def get_eta_s():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function

    """cello in gd
    @legacy_function      
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function
    """
    """cello in gd
    @legacy_function      
    def get_potential_energy():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function
    """

    @legacy_function      
    def get_energy_error():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function

    def get_energies(self):
        energy_unit = nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2
        kinetic_energy = self.get_kinetic_energy() | energy_unit
        potential_energy = self.get_potential_energy() | energy_unit
        return (self.convert_nbody.to_si(kinetic_energy), self.convert_nbody.to_si(potential_energy))

    @legacy_function      
    def find_colliding_secondary():
        function = LegacyFunctionSpecification()  
        function.addParameter('id1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    """
    @legacy_function          
    def remove_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function    
    """    
        
class PhiGRAPEBinding(NBodyGravitationalDynamicsBinding):
   
    parameter_definitions = [
        parameters.ModuleMethodParameterDefinition(
            "get_eps2", "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        ),
        parameters.ModuleMethodParameterDefinition(
            "get_eta", "set_eta1",
            "eta", 
            "timestep parameter", 
            units.none , 
            0.01 |  units.none
        ),
        parameters.ModuleMethodParameterDefinition(
            "get_eta_s", "set_eta_s",
            "eta_s", 
            "parameter to determine the initial timestep", 
            units.none , 
            0.002 |  units.none
        ),
    ]

    attribute_definitions = [
        attributes.AttributeDefinition(
            name = "mass",
            setup_parameters = ["mass"],
            setter = ("set_mass", ["mass"]),
            description = "mass of a star",
            unit = nbody_system.mass,
            default = 1 | nbody_system.mass          
        ),
        attributes.AttributeDefinition(
            name = "radius",
            setup_parameters = ["radius"],
            setter = ("set_radius", ["radius"]),
            description = "radius of a star",
            unit = nbody_system.length,
            default = 1 | nbody_system.length          
        ),
        attributes.AttributeDefinition(
            names = ["x","y","z"],
            setup_parameters = ["x","y","z"],
            setter = ("set_position", ["x","y","z"]),
            description = "coordinate of a star",
            unit = nbody_system.length,
            default = 0.0 | nbody_system.length          
        ),
        attributes.AttributeDefinition(
            names = ["vx","vy","vz"],
            setup_parameters = ["vx","vy","vz"],
            setter = ("set_velocity", ["vx","vy","vz"]),
            description = "coordinate of a star",
            unit = nbody_system.speed,
            default = 0.0 | nbody_system.speed          
        ),
    ]


    def __init__(self, convert_nbody = None):
        NBodyGravitationalDynamicsBinding.__init__(self, convert_nbody)
    
    def current_model_time(self):
        return self.convert_nbody.to_si( self.t | nbody_system.time)
            

class PhiGRAPE(PhiGRAPEInterface, PhiGRAPEBinding):
    
    def __init__(self, convert_nbody = None):
        PhiGRAPEInterface.__init__(self)
        PhiGRAPEBinding.__init__(self, convert_nbody)
        
  
