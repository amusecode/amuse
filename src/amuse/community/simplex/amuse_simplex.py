from amuse.community import *

class SimpleX(LegacyInterface):
    include_headers=['worker.h']
    
    def __init__(self, **kwargs):       
        LegacyInterface.__init__(self,name_of_the_worker="simplex_worker",**kwargs)
        

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
    def initialize():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        function.addParameter('current_time', dtype='d', direction=function.IN)
        return function;

    @legacy_function    
    def reinitialize():
        function = LegacyFunctionSpecification() 
        function.result_type = 'i'
        return function;
        
    @legacy_function    
    def add_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['x','y','z','rho','flux','xion']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z','rho','flux','xion']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = None
        return function
        
    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.addParameter('synchronize', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def remove_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
