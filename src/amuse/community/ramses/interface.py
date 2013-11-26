from amuse.community import *
from amuse.community.interface.hydro import HydrodynamicsInterface
from amuse.community.interface.common import CommonCode
from amuse.units.generic_unit_system import *

class RamsesInterface(CodeInterface, HydrodynamicsInterface, 
        StoppingConditionInterface, LiteratureReferencesMixIn, CodeWithDataDirectories):
    """
    Ramses is an AMR hydrodynamics code...
    
    The reference for Ramses (original version):
        .. [#] ...
    """
    
    MODE_3D = '3d'
    MODE_2D = '2d'
    MODE_1D = '1d'
    modes = [MODE_3D, MODE_2D, MODE_1D]
    
    def __init__(self, number_of_workers=1, mode=MODE_3D, **options):
        CodeInterface.__init__(self, 
            name_of_the_worker=self.name_of_the_worker(mode), 
            number_of_workers=number_of_workers, **options)
        LiteratureReferencesMixIn.__init__(self)
    
    def name_of_the_worker(self, mode):
        if mode == self.MODE_3D:
            return 'ramses_worker'
        elif mode == self.MODE_2D:
            return 'ramses_worker_2d'
        elif mode == self.MODE_1D:
            return 'ramses_worker_1d'
        else:
            return 'ramses_worker'
    
#    @legacy_function
#    def set_ramses_data_directory():
#        function = LegacyFunctionSpecification()  
#        function.addParameter('data_directory', dtype='string', direction=function.IN,
#            description = "Name of the data directory")
#        function.result_type = 'int32'
#        return function
    
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
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('number_of_points', 'i', function.LENGTH)            
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_grid_state():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
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
    def evolve_model():
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
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_index_of_position():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
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
    def get_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx','rhovy','rhovz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    

    @legacy_function
    def get_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['energy']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    
    def get_number_of_grids():
        return (0,1)
    
    
    set_grid_energy_density = None
    set_grid_density = None
    set_grid_momentum_density = None
    
    
    
    @legacy_function
    def get_boundary_index_range_inclusive():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN)
        function.addParameter('minx', dtype='i', direction=function.OUT)
        function.addParameter('maxx', dtype='i', direction=function.OUT)
        function.addParameter('miny', dtype='i', direction=function.OUT)
        function.addParameter('maxy', dtype='i', direction=function.OUT)
        function.addParameter('minz', dtype='i', direction=function.OUT)
        function.addParameter('maxz', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_boundary_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_boundary_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN, default = 1)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    
    @legacy_function    
    def get_boundary_position_of_index():
        """
        Retrieves the x, y and z position of the center of
        the cell with coordinates i, j, k 
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN, default = 1)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_parallel_decomposition():
        """Set the number of processes per dimension, the
        total number of available processors must be
        nx * ny * nz
        """
        function = LegacyFunctionSpecification()
        for x in ['nx','ny','nz']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_parallel_decomposition():
        """Retrieves the number of processes per dimension, the
        total number of available processors must be
        nx * ny * nz
        """
        function = LegacyFunctionSpecification()
        for x in ['nx','ny','nz']:
            function.addParameter(x, dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_timestep():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.OUT) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_timestep():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.IN) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_gamma():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_gamma():
        function = LegacyFunctionSpecification()  
        function.addParameter('value', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function    
    def get_hydro_state_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, default = 0)
        for x in ['rho','rhovx','rhovy','rhovz','rhoe']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'i' 
        function.must_handle_array = True
        return function
    

class Ramses(CommonCode):

    modes = RamsesInterface.modes
        
    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        self.stopping_conditions = StoppingConditions(self)
        
        CommonCode.__init__(self,  RamsesInterface(**options), **options)
#        self.set_ramses_data_directory(self.data_directory)
    

    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
            
    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")
        
    def define_methods(self, object):
        object.add_method('evolve_model', (time,), (object.ERROR_CODE,))
        object.add_method(
            'get_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX,),
            (length, length, length, object.ERROR_CODE,)
        )
        
        density = mass / (length**3)
        momentum_density =  mass / (time * (length**2))
        energy_density =  mass / ((time**2) * length)
        acceleration =  length / time ** 2
        
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
            'get_grid_density',
            (object.INDEX, object.INDEX, object.INDEX),
            (density,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX),
            (momentum_density, momentum_density, momentum_density,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_energy_density',
            (object.INDEX, object.INDEX, object.INDEX),
            (energy_density,
            object.ERROR_CODE,)
        )
    
        object.add_method(
            'get_gravity_field',
            (object.INDEX, object.INDEX, object.INDEX),
            (acceleration, acceleration, acceleration,
            object.ERROR_CODE,)
        )
        object.add_method(
            'set_gravity_field',
            (object.INDEX, object.INDEX, object.INDEX,
             acceleration, acceleration, acceleration,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_time',
            (),
            (time, object.ERROR_CODE,)
        )
    
        object.add_method(
            'setup_mesh',
            (object.NO_UNIT, object.NO_UNIT, object.NO_UNIT, length, length, length,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'set_boundary',
            (object.NO_UNIT, object.NO_UNIT, object.NO_UNIT, object.NO_UNIT, object.NO_UNIT, object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'set_parallel_decomposition',
            (object.NO_UNIT, object.NO_UNIT, object.NO_UNIT),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_boundary_state',
            (object.INDEX, object.INDEX, object.INDEX,
            density, momentum_density, momentum_density, momentum_density, energy_density,
            object.INDEX),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_boundary_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density, momentum_density, momentum_density, momentum_density, energy_density,
            object.ERROR_CODE,)
        )
        object.add_method(
            'get_boundary_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length, length, length, object.ERROR_CODE,)
        )
        object.add_method(
            'get_boundary_index_range_inclusive',
            (object.INDEX),
            (object.NO_UNIT, object.NO_UNIT,object.NO_UNIT, object.NO_UNIT,object.NO_UNIT, object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_timestep",
            (),
            (time, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_timestep",
            (time, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_hydro_state_at_point',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,
                generic_unit_system.speed, generic_unit_system.speed, generic_unit_system.speed),
            (generic_unit_system.density, generic_unit_system.momentum_density, generic_unit_system.momentum_density, 
                generic_unit_system.momentum_density, generic_unit_system.energy_density, object.ERROR_CODE)
        )
    
        self.stopping_conditions.define_methods(object)
    
    def define_particle_sets(self, object):
        object.define_grid('grid')
        object.set_grid_range('grid', 'get_index_range_inclusive')
        object.add_getter('grid', 'get_position_of_index', names=('x','y','z'))
    
        object.add_getter('grid', 'get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        object.add_setter('grid', 'set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        
        object.add_getter('grid', 'get_grid_density', names=('rho',))
        object.add_getter('grid', 'get_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        object.add_getter('grid', 'get_grid_energy_density', names=('energy',))
        
        object.define_grid('acceleration_grid')
        object.set_grid_range('acceleration_grid', 'get_index_range_inclusive')
        object.add_getter('acceleration_grid', 'get_position_of_index', names=('x','y','z'))
        object.add_setter('acceleration_grid', 'set_gravity_field', names=('ax','ay','az'))
        object.add_getter('acceleration_grid', 'get_gravity_field', names=('ax','ay','az'))
    
        #object.add_setter('grid', 'set_momentum_density', names=('rhovx','rhovy','rhovz'))
        #object.add_setter('grid', 'set_density', names=('rho',))
        #object.add_setter('grid', 'set_energy_density', names=('energy',))
        
        
        

    def define_parameters(self, object):
        
        
        object.add_method_parameter(
            "get_gamma", 
            "set_gamma",
            "gamma", 
            "ratio of specific heats used in equation of state", 
            default_value = 1.6666666666666667
        )
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshx",
            "nx", 
            "number of cells in the x direction",
            10,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshy",
            "ny", 
            "number of cells in the y direction",
            10,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshz",
            "nz", 
            "number of cells in the z direction",
            10,
        )
        
        object.add_caching_parameter(
            "setup_mesh", 
            "xlength",
            "length_x", 
            "length of model in the x direction", 
            10 | length,
        )
        object.add_caching_parameter(
            "setup_mesh", 
            "ylength",
            "length_y", 
            "length of model in the x direction",
            10 | length,
        )
        object.add_caching_parameter(
            "setup_mesh", 
            "zlength",
            "length_z", 
            "length of model in the z direction",
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
            "set_parallel_decomposition", 
            "nx",
            "nproc_x", 
            "number of processors for the x direction",
            0,
        )
        object.add_caching_parameter(
            "set_parallel_decomposition", 
            "ny",
            "nproc_y", 
            "number of processors for the y direction",
            0,
        )
        object.add_caching_parameter(
            "set_parallel_decomposition", 
            "nz",
            "nproc_z", 
            "number of processors for the z direction",
            0,
        )
        
        object.add_vector_parameter(
            "parallel_decomposition",
            "number of processors for each dimensions",
            ("nproc_x", "nproc_y", "nproc_z")
        )
        
        object.add_caching_parameter(
             "set_boundary", 
             "xbound1",
             "xbound1", 
             "boundary conditions on first (inner, left) X boundary", 
             "reflective",
        )
        object.add_caching_parameter(
            "set_boundary", 
            "xbound2",
            "xbound2", 
            "boundary conditions on second (outer, right) X boundary",
            "reflective",
        )
        
        object.add_caching_parameter(
            "set_boundary", 
            "ybound1",
            "ybound1", 
            "boundary conditions on first (inner, front) Y boundary",
            "reflective",
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "ybound2",
            "ybound2", 
            "boundary conditions on second (outer, back) Y boundary",
            "reflective",
        )
        
        object.add_caching_parameter(
            "set_boundary", 
            "zbound1",
            "zbound1", 
            "boundary conditions on first (inner, bottom) Z boundary",
            "reflective",
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "zbound2",
            "zbound2", 
            "boundary conditions on second (outer, top) Z boundary",
            "reflective",
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
        
        self.stopping_conditions.define_parameters(object)
    
    def get_index_range_inclusive(self):
        """
        Returns the min and max values of indices in each
        direction. The range is inclusive, the min index
        and max index both exist and can be queried.
        The total number of cells in one direction
        is max - min + 1.
        """
        nx, ny, nz = self.get_mesh_size()
        return (1, nx, 1, ny, 1, nz)
    

    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
        
    def itergrids(self):
        yield self.grid
    
    
    
    def define_state(self, object): 
        CommonCode.define_state(self, object)   
        #object.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        
        object.add_method('INITIALIZED', 'set_parallel_decomposition')
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        
        object.add_transition('EDIT', 'RUN', 'initialize_grid')
        object.add_method('RUN', 'evolve_model')
        object.add_method('RUN', 'get_hydro_state_at_point')
        
        for state in ['EDIT', 'RUN']:
            for methodname in [
                    'get_grid_state',
                    'set_grid_state',
                    'get_grid_density',
                    'set_grid_density',
                    'set_grid_energy_density',
                    'get_grid_energy_density',
                    'get_grid_momentum_density',
                    'set_grid_momentum_density', 
                    'get_position_of_index',
                    'get_index_of_position',
                    'set_grid_scalar',
                    'get_grid_scalar',
                    'get_mesh_size',
                    'set_gravity_field',
                    'get_gravity_field',
                    'get_boundary_state',
                    'set_boundary_state',
                    'get_boundary_position_if_index',
                    'get_boundary_index_range_inclusive'
                ]:
                object.add_method(state, methodname)    
    
    def specify_boundary_grid(self, definition, index_of_boundary, index_of_grid = 1):
        definition.set_grid_range('get_boundary_index_range_inclusive')
        
        definition.add_getter('get_boundary_position_of_index', names=('x','y','z'))
        
        definition.add_getter('get_boundary_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_boundary_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
       
        definition.define_extra_keywords({'index_of_boundary': index_of_boundary})
    
    BOUNDARY_NAME_TO_INDEX = {
        'xbound1': 1,
        'xbound2': 2,
        'ybound1': 3,
        'ybound2': 4,
        'zbound1': 5,
        'zbound2': 6,
    }
    def get_boundary_grid(self, name):
        if not name in self.BOUNDARY_NAME_TO_INDEX:
            raise Exception("boundary name is not known {0}".format(name))
        index_of_boundary = self.BOUNDARY_NAME_TO_INDEX[name]
        
        return self._create_new_grid(self.specify_boundary_grid, index_of_boundary = index_of_boundary, index_of_grid = 1)
    
