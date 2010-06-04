from amuse.legacy import *
from amuse.legacy.interface.common import CommonCodeInterface

class AthenaInterface(LegacyInterface, CommonCodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        LegacyInterface.__init__(self, **options)
    
    @legacy_function
    def test():
        function = LegacyFunctionSpecification()  
        #function.addParameter('mass', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
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
    def par_seti():
        function = LegacyFunctionSpecification() 
        function.addParameter('block', dtype='s', direction=function.IN) 
        function.addParameter('name', dtype='s', direction=function.IN) 
        function.addParameter('fmt', dtype='s', direction=function.IN) 
        function.addParameter('ival', dtype='int32', direction=function.IN) 
        function.addParameter('comment', dtype='s', direction=function.IN) 
        function.result_type = None
        return function
        
    @legacy_function
    def par_geti():
        function = LegacyFunctionSpecification() 
        function.addParameter('block', dtype='s', direction=function.IN) 
        function.addParameter('name', dtype='s', direction=function.IN) 
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def par_setd():
        function = LegacyFunctionSpecification() 
        function.addParameter('block', dtype='s', direction=function.IN) 
        function.addParameter('name', dtype='s', direction=function.IN) 
        function.addParameter('fmt', dtype='s', direction=function.IN) 
        function.addParameter('dval', dtype='float64', direction=function.IN) 
        function.addParameter('comment', dtype='s', direction=function.IN) 
        function.result_type = None
        return function
        
    @legacy_function
    def par_getd():
        function = LegacyFunctionSpecification() 
        function.addParameter('block', dtype='s', direction=function.IN) 
        function.addParameter('name', dtype='s', direction=function.IN) 
        function.result_type = 'float64'
        return function
          
    def setup_mesh(self, nmeshx, nmeshy, nmeshz, xlength, ylength, zlength):
        self.par_seti("grid", "Nx1", "%d", nmeshx, "")
        self.par_seti("grid", "Nx2", "%d", nmeshy, "")
        self.par_seti("grid", "Nx3", "%d", nmeshz, "")
        
        self.par_setd("grid", "x1min", "%.15e", 0.0, "")
        self.par_setd("grid", "x1max", "%.15e", xlength, "")
        self.par_setd("grid", "x2min", "%.15e", 0.0, "")
        self.par_setd("grid", "x2max", "%.15e", ylength, "")
        self.par_setd("grid", "x3min", "%.15e", 0.0, "")
        self.par_setd("grid", "x3max", "%.15e", zlength, "")
        
        return 0
        
    
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
        
    def set_isocsound(self, value):
        self.par_setd("problem", "iso_csound", "%.15e", value, "") 
        
    def set_gamma(self, value):
        self.par_setd("problem", "gamma", "%.15e", value, "") 
    
    def set_courant_friedrichs_lewy_number(self, value):
        self.par_setd("time", "cour_no", "%.15e", value, "") 
        
    def set_boundary(self, xbound1, xbound2, ybound1, ybound2, zbound1, zbound2):
        map_from_string_to_flag = {"reflective": 1, "outflow":2, "periodic":4}
        
        self.par_seti("grid", "ibc_x1", "%d", map_from_string_to_flag[xbound1], "")
        self.par_seti("grid", "obc_x1", "%d", map_from_string_to_flag[xbound2], "")
        self.par_seti("grid", "ibc_x2", "%d", map_from_string_to_flag[ybound1], "")
        self.par_seti("grid", "obc_x2", "%d", map_from_string_to_flag[ybound2], "")
        self.par_seti("grid", "ibc_x3", "%d", map_from_string_to_flag[zbound1], "")
        self.par_seti("grid", "obc_x3", "%d", map_from_string_to_flag[zbound1], "")
        
    
    @legacy_function    
    def initialize_grid():
        function = LegacyFunctionSpecification()  
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
    def get_nghost():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='int32', direction=function.OUT) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.OUT) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def esys_roe_adb_hydro():
        function = LegacyFunctionSpecification() 
        function.addParameter('index', dtype='int32', direction=function.OUT) 
        function.addParameter('u', dtype='float64', direction=function.OUT) 
        function.addParameter('v', dtype='float64', direction=function.OUT) 
        function.addParameter('w', dtype='float64', direction=function.OUT) 
        function.addParameter('h', dtype='float64', direction=function.OUT) 
        function.addParameter('ev', dtype='float64', direction=function.OUT) 
        for i in range(5):
            function.addParameter('rem{0}'.format(i), dtype='float64', direction=function.OUT) 
        for i in range(5):
            function.addParameter('lem{0}'.format(i), dtype='float64', direction=function.OUT) 
        function.result_type = 'i'
        function.can_handle_array = True
        return function
        
    @legacy_function
    def fill_grid_linearwave_1d():
        function = LegacyFunctionSpecification() 
        function.addParameter('wave_flag', dtype='int32', direction=function.IN) 
        function.addParameter('amplitude', dtype='float64', direction=function.IN) 
        function.addParameter('vflow', dtype='float64', direction=function.IN) 
        function.addParameter('wave_dir', dtype='int32', direction=function.IN) 
        function.result_type = 'i'
        function.can_handle_array = True
        return function
        
    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='float64', direction=function.IN) 
        function.result_type = 'i'
        return function
        
    
class Athena(CodeInterface):

    def __init__(self, **options):
        CodeInterface.__init__(self,  AthenaInterface(**options), **options)
    
    def define_methods(self, object):
        pass
    
    
    def define_particle_sets(self, object):
        pass
        
    def _evolve_particles(self, particles, end_time):
        pass

    def evolve_model(self, end_time = None):
        pass
