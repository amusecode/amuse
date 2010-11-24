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
    
    
    def __init__(self, number_of_workers = 1, **options):
        LegacyInterface.__init__(self, self.name_of_the_worker(number_of_workers),**options)
        LiteratureRefs.__init__(self)
    
    def name_of_the_worker(self, number_of_workers):
        if number_of_workers > 1:
            return 'capreole_worker_mpi'
        else:
            return 'capreole_worker'
    
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
        #function.addParameter('time', dtype='d', direction=function.IN)
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
    def set_grid_state():
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
    
    
    

    

    def get_index_range_inclusive(self):
        ni,nj,nk,error = self.get_mesh_size()
        return (1, ni, 1, nj, 1, nk)
    
    

    

    

    @legacy_function
    def get_momentum_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx','rhovy','rhovz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    

    

    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    

    @legacy_function
    def get_energy_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['energy']:
            function.addParameter(x, dtype='d', direction=function.OUT)
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
            (object.INDEX, object.INDEX, object.INDEX,),
            (length, length, length, object.ERROR_CODE,)
        )
        
        density = mass / (length**3)
        momentum_density =  mass / (time * (length**2))
        energy_density =  mass / ((time**2) * length)
        
        object.add_method(
            'set_grid_state',
            (object.INDEX, object.INDEX, object.INDEX,
            density, momentum_density, momentum_density, momentum_density, energy_density,
            ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_grid_state',
            (object.INDEX, object.INDEX, object.INDEX),
            (density, momentum_density, momentum_density, momentum_density, energy_density,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_density',
            (object.INDEX, object.INDEX, object.INDEX),
            (density,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX),
            (momentum_density, momentum_density, momentum_density,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_energy_density',
            (object.INDEX, object.INDEX, object.INDEX),
            (energy_density,
            object.ERROR_CODE,)
        )
    
        
    
    def define_particle_sets(self, object):
        object.define_grid('grid')
        object.set_grid_range('grid', 'get_index_range_inclusive')
        object.add_getter('grid', 'get_position_of_index', names=('x','y','z'))
    
        object.add_getter('grid', 'get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        object.add_setter('grid', 'set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        
        object.add_getter('grid', 'get_density', names=('rho',))
        object.add_getter('grid', 'get_momentum_density', names=('rhovx','rhovy','rhovz'))
        object.add_getter('grid', 'get_energy_density', names=('energy',))
        
        #object.add_setter('grid', 'set_momentum_density', names=('rhovx','rhovy','rhovz'))
        #object.add_setter('grid', 'set_density', names=('rho',))
        #object.add_setter('grid', 'set_energy_density', names=('energy',))
        
        
        

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
    
        object.add_caching_parameter(
            "set_boundary", 
            "xbound1",
            "xbound1", 
            "boundary conditions on first (inner, left) X boundary", 
            units.string, 
            "reflective" | units.string,
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "xbound2",
            "xbound2", 
            "boundary conditions on second (outer, right) X boundary", 
            units.string, 
            "reflective" | units.string,
        )
        
        object.add_caching_parameter(
            "set_boundary", 
            "ybound1",
            "ybound1", 
            "boundary conditions on first (inner, front) Y boundary", 
            units.string, 
            "reflective" | units.string,
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "ybound2",
            "ybound2", 
            "boundary conditions on second (outer, back) Y boundary", 
            units.string, 
            "reflective" | units.string,
        )
        
        object.add_caching_parameter(
            "set_boundary", 
            "zbound1",
            "zbound1", 
            "boundary conditions on first (inner, bottom) Z boundary", 
            units.string, 
            "reflective" | units.string,
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "zbound2",
            "zbound2", 
            "boundary conditions on second (outer, top) Z boundary", 
            units.string, 
            "reflective" | units.string,
        )
        
        
        
        object.add_vector_parameter(
            "x_boundary_conditions",
            "boundary conditions for the X directorion",
            ("xbound1", "xbound2")
        )
        
        
        object.add_vector_parameter(
            "y_boundary_conditions",
            "boundary conditions for the Y directorion",
            ("ybound1", "ybound2")
        )
        
        
        object.add_vector_parameter(
            "z_boundary_conditions",
            "boundary conditions for the Z directorion",
            ("zbound1", "zbound2")
        )


    def get_index_range_inclusive(self):
        """
        Returns the min and max values of indices in each
        direction. The range is inclusive, the min index
        and max index both exist and can be queried.
        The total number of cells in one direction
        is max - min + 1.
        """
        nx, ny, nz, error = self.get_mesh_size()
        return (1, nx, 1, ny, 1, nz)
    
    

    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
    
    
