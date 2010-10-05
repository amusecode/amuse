from amuse.legacy import *

from amuse.support.units.generic_unit_system import *
from amuse.legacy.interface.common import CommonCodeInterface

class CapreoleInterface(LegacyInterface, CommonCodeInterface, LiteratureRefs):
    """
    Capreole is a grid-based astrophysical hydrodynamics code developed by Garrelt Mellema. 
    It works in one, two dimensions, and three spatial dimensions and is programmed in 
    Fortran 90. It is parallelized with MPI. For the hydrodynamics it relies on the 
    Roe-Eulderink-Mellema (REM) solver, which is an approximate Riemann solver for arbitrary
    metrics. It can solve different hydrodynamics problems. Capreole has run on single 
    processors, but also on massively parallel systems (e.g. 512 processors on a BlueGene/L).
    
    The reference for Capreole (original version):
        .. [#] Mellema, Eulderink & Icke 1991, A&A 252, 718
    """
    
    def __init__(self,name_of_the_worker = 'worker',**args):
        LegacyInterface.__init__(self, name_of_the_worker,**args)
        LiteratureRefs.__init__(self)
    
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
        
        


    @legacy_function
    def get_mesh_size():
        function = LegacyFunctionSpecification()
        function.addParameter('nmeshx', dtype='i', direction=function.OUT)
        function.addParameter('nmeshy', dtype='i', direction=function.OUT)
        function.addParameter('nmeshz', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    
    

    
class GLCapreoleInterface(CapreoleInterface):
    def __init__(self, **options):
        LegacyInterface.__init__(self,name_of_the_worker = 'glworker', **options)
        
    @legacy_function
    def viewer():
        function = LegacyFunctionSpecification()
        return function


class Capreole(CodeInterface):

    def __init__(self, **options):
        CodeInterface.__init__(self,  CapreoleInterface(**options), **options)
    
    def define_properties(self, object):
        object.add_property('get_time', time, "model_time")
        
    def define_methods(self, object):
        object.add_method(
            'evolve',
            (time,),
            (object.ERROR_CODE,),
        )
        object.add_method(
            'get_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length, length, length, object.ERROR_CODE,)
        )
        
        density = mass / (length**3)
        momentum =  mass / (time * (length**2))
        energy =  mass / ((time**2) * length)
        
        object.add_method(
            'fill_grid_state',
            (object.INDEX, object.INDEX, object.INDEX,
            density, momentum, momentum, momentum, energy,
            object.INDEX, object.INDEX),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_grid_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density, momentum, momentum, momentum, energy,
            object.ERROR_CODE,)
        )
    
        
    
    def define_particle_sets(self, object):
        object.define_grid('grid')
        object.set_grid_range('grid', 'get_index_range_inclusive')
        object.add_getter('grid', 'get_position_of_index', names=('x','y','z'))
        object.add_getter('grid', 'get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        object.add_setter('grid', 'fill_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
    
        
        

    def define_parameters(self, object):
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshx",
            "nx", 
            "number of cells in the x direction", 
            units.none, 
            10 | units.none,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshy",
            "ny", 
            "number of cells in the y direction", 
            units.none, 
            10 | units.none,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshz",
            "nz", 
            "number of cells in the z direction", 
            units.none, 
            10 | units.none,
        )
        
        object.add_caching_parameter(
            "setup_mesh", 
            "xlength",
            "length_x", 
            "length of model in the x direction", 
            length, 
            10 | length,
        )
        object.add_caching_parameter(
            "setup_mesh", 
            "ylength",
            "length_y", 
            "length of model in the x direction", 
            length, 
            10 | length,
        )
        object.add_caching_parameter(
            "setup_mesh", 
            "zlength",
            "length_z", 
            "length of model in the z direction", 
            length, 
            10 | length,
        )
        
        object.add_vector_parameter(
            "mesh_size",
            "number of cells in the x, y and z directions",
            ("nx", "ny", "nz")
        )
        
        object.add_vector_parameter(
            "mesh_length",
            "length of the model in the x, y and z directions",
            ("length_x", "length_y", "length_z")
        )


    def get_index_range_inclusive(self):
        """
        Returns the min and max values of indices in each
        direction. The range is inclusive, the min index
        and max index both exist and can be queried.
        The total number of cells in one direction
        is max - min + 1.
        """
        nx, ny, nz = self.get_mesh_Size()
        return (1, nx, 1, ny, 1, nz)
    
    
