import numpy

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.legacy import *
from amuse.legacy.interface.gd import GravitationalDynamics
from amuse.legacy.interface.gd import GravitationalDynamicsInterface
from amuse.legacy.support.lit import LiteratureRefs
from amuse.support.data import binding

class PhiGRAPEInterface(LegacyInterface, LiteratureRefs, GravitationalDynamics):
    """
        .. [#] Refs not included yet
    """    
    
    MODE_G6LIB = 'g6lib'
    MODE_GPU   = 'gpu'
    MODE_GRAPE = 'grape'
    MODE_PG    = 'pg'
    
    def __init__(self, convert_nbody = None, mode = MODE_G6LIB):
        LegacyInterface.__init__(self, name_of_the_worker = self.name_of_the_muse_worker(mode))
        LiteratureRefs.__init__(self)

    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'worker_code'
        elif mode == self.MODE_GPU:
            return 'worker_code_gpu'
        elif mode == self.MODE_GRAPE:
            return 'worker_code_grape'
        elif mode == self.MODE_PG:
            return 'worker_code_phantom_grape'
            
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
        function.result_type = 'int32'
        return function

    @legacy_function      
    def set_eta1():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type = 'int32'
        return function


    @legacy_function      
    def get_eta():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='d', direction=function.OUT) 
        function.result_type = 'int32'
        return function
        
    @legacy_function      
    def get_eta_s():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type = 'int32'
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


            

class PhiGRAPEInterfaceGL(PhiGRAPEInterface):
    
    def __init__(self, mode = PhiGRAPEInterface.MODE_G6LIB):
        PhiGRAPEInterface.__init__(self, mode = mode)

    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()  
        return function
  
    def name_of_the_muse_worker(self, mode):
        if mode == self.MODE_G6LIB:
            return 'glworker_code'
        if mode == self.MODE_GPU:
            return 'glworker_code_gpu'
        if mode == self.MODE_GRAPE:
            return 'glworker_code_grape'
            
            

class PhiGRAPE(GravitationalDynamicsInterface):
    
    
    def __init__(self, convert_nbody = None, mode = PhiGRAPEInterface.MODE_G6LIB, use_gl = False):
        
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
       
        nbody_interface = None
        if use_gl:
            nbody_interface = PhiGRAPEInterfaceGL(mode)
        else:
            nbody_interface = PhiGRAPEInterface(mode)
        
        GravitationalDynamicsInterface.__init__(
            self,
            nbody_interface,
            convert_nbody,
        )     
            
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2", 
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        ) 
        
        object.add_method_parameter(
            "get_eta", 
            "set_eta1",
            "timestep_parameter", 
            "timestep parameter", 
            units.none , 
            0.02 |  units.none
        )
        
        object.add_method_parameter(
            "get_eta_s", 
            "set_eta_s",
            "initial_timestep_parameter", 
            "parameter to determine the initial timestep", 
            units.none , 
            0.01 |  units.none
        )
  
           



