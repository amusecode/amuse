from amuse.legacy import *


class Capreole(LegacyInterface):
    def __init__(self,name_of_the_worker = 'worker',**args):
        LegacyInterface.__init__(self, name_of_the_worker,**args)

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
    def setup_mesh():
        function = LegacyFunctionSpecification()  
        function.addParameter('nmeshx', dtype='i', direction=function.IN)
        function.addParameter('nmeshy', dtype='i', direction=function.IN)
        function.addParameter('nmeshz', dtype='i', direction=function.IN)
        function.addParameter('xlength', dtype='d', direction=function.IN)
        function.addParameter('ylength', dtype='d', direction=function.IN)
        function.addParameter('zlength', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function


    @legacy_function    
    def initialize_grid():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary():
        function = LegacyFunctionSpecification()  
        for x in ["xbound1","xbound2","ybound1","ybound2","zbound1","zbound2"]:
            function.addParameter(x, dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary_innerxstate():
        function = LegacyFunctionSpecification()  
        for x in ["rho_in","rhvx_in","rhvy_in","rhvz_in","en_in"]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary_outerxstate():
        function = LegacyFunctionSpecification()  
        for x in ["rho_in","rhvx_in","rhvy_in","rhvz_in","en_in"]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary_innerystate():
        function = LegacyFunctionSpecification()  
        for x in ["rho_in","rhvx_in","rhvy_in","rhvz_in","en_in"]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary_outerystate():
        function = LegacyFunctionSpecification()  
        for x in ["rho_in","rhvx_in","rhvy_in","rhvz_in","en_in"]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary_innerzstate():
        function = LegacyFunctionSpecification()  
        for x in ["rho_in","rhvx_in","rhvy_in","rhvz_in","en_in"]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_boundary_outerzstate():
        function = LegacyFunctionSpecification()  
        for x in ["rho_in","rhvx_in","rhvy_in","rhvz_in","en_in"]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def fill_grid_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_grid_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_gravity_field():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['f_x','f_y','f_z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_gravity_field():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['f_x','f_y','f_z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function    
    def get_position_of_index():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['i','j','k']:
          function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['x','y','z']:
          function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_index_of_position():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['x','y','z']:
          function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['i','j','k']:
          function.addParameter(x, dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function


class GLCapreole(Capreole):
    def __init__(self, **options):
        LegacyInterface.__init__(self,name_of_the_worker = 'glworker', **options)
        
    @legacy_function
    def viewer():
        function = LegacyFunctionSpecification()
        return function 
