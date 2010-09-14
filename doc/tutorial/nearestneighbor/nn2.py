from amuse.legacy import *

class NearestNeighborInterface(LegacyInterface):
    
    use_modules = ['NN']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, **keyword_arguments)
    
    @legacy_function
    def commit_parameters():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_maximum_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_maximum_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function    

    @legacy_function
    def run():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function    
    
    @legacy_function
    def get_close_neighbors():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('index_of_first_neighbor', dtype='float64', direction=function.OUT)
        function.addParameter('index_of_second_neighbor', dtype='float64', direction=function.OUT)
        function.addParameter('index_of_third_neighbor', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    
    @legacy_function
    def get_nearest_neighbor():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_neighbor', dtype='float64', direction=function.OUT)
        function.addParameter('distance', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    
    
class NearestNeighbor(CodeInterface):

    def __init__(self):
        CodeInterface.__init__(self,  NearestNeighborInterface())
    
