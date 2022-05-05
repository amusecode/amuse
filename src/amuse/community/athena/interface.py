from amuse.community import *

from amuse.community.interface.mhd import MagnetohydrodynamicsInterface
from amuse.community.interface.common import CommonCode

import numpy

from amuse.units.generic_unit_system import *

class AthenaInterface(CodeInterface, MagnetohydrodynamicsInterface, LiteratureReferencesMixIn, StoppingConditionInterface):
    """
    Athena is a grid-based code for astrophysical hydrodynamics. Athena can solve 
    magnetohydrodynamics (MHD) as well, but this is currently not supported from 
    AMUSE. It was developed primarily for studies of the interstellar medium, 
    star formation, and accretion flows.
    
    The current version (Athena v4.0) implements algorithms for the following physics:
    * compressible hydrodynamics and MHD in 1D, 2D, and 3D,
    * ideal gas equation of state with arbitrary gamma (including gamma = 1, an isothermal EOS),
    * an arbitrary number of passive scalars advected with the flow,
    * self-gravity, and/or a static gravitational potential,
    * Ohmic resistivity, ambipolar diffusion, and the Hall effect,
    * both Navier-Stokes and anisotropic (Braginskii) viscosity,
    * both isotropic and anisotropic thermal conduction,
    * optically-thin radiative cooling. 
    
    In addition, Athena allows for the following grid and parallelization options:
    * Cartesian or cylindrical coordinates,
    * static (fixed) mesh refinement,
    * shearing-box source terms, and an orbital advection algorithm for MHD,
    * parallelization using domain decomposition and  MPI. 
    
    A variety of choices are also available for the numerical algorithms, such as 
    different Riemann solvers and spatial reconstruction methods.
    
    The relevant references are:
        .. [#] Gardiner & Stone 2005, JCP, 205, 509  (2D JCP Method) [2005JCoPh.205..509G]
        .. [#] Gardiner & Stone 2007, JCP, 227, 4123 (3D JCP Method) [2008JCoPh.227.4123G]
        .. [#] Stone et al. 2008, ApJS, 178, 137 (Method) [2008ApJS..178..137S]
        .. [#] Stone & Gardiner 2009, NewA, 14, 139 (van Leer Integrator) [2009NewA...14..139S]
        .. [#] Skinner & Ostriker 2010, ApJ, 188, 290 (Cylindrical Integrator) [2010ApJS..188..290S]
        .. [#] Stone & Gardiner 2010, ApJS, 189, 142 (Shearing Box Method) [2010ApJS..189..142S]
    """
    
    include_headers = ['worker_code.h', 'stopcond.h']
    
    MODE_NORMAL = 'normal'
    MODE_SELF_GRAVITY   = 'self-gravity'
    MODE_MHD   = 'mhd'
    MODE_SCALAR   = 'scalar'
    
    def __init__(self, mode = MODE_NORMAL, **options):
        
        self.mode = mode
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **options)
        self.set_auto_decomposition(1)
        self.par_seti("domain1", "AutoWithNProc", "%d", self.channel.number_of_workers, "-")
        LiteratureReferencesMixIn.__init__(self)
        self.number_of_grids = 1
        
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'athena_worker'
        elif mode == self.MODE_SELF_GRAVITY:
            return 'athena_worker_selfgrav'
        elif mode == self.MODE_MHD:
            return 'athena_worker_mhd'
        elif mode == self.MODE_SCALAR:
            return 'athena_worker_scalar'
        else:
            return 'athena_worker'
        
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
        self.par_seti("job","num_domains", "%d", self.number_of_grids, "-")
        self.par_seti("domain1", "level", "%d", 0, "-")
        
        self.par_seti("domain1", "Nx1", "%d", nmeshx, "-")
        self.par_seti("domain1", "Nx2", "%d", nmeshy, "-")
        self.par_seti("domain1", "Nx3", "%d", nmeshz, "-")
        
        self.par_setd("domain1", "x1min", "%.15e", 0.0, "-")
        self.par_setd("domain1", "x1max", "%.15e", xlength, "-")
        self.par_setd("domain1", "x2min", "%.15e", 0.0, "-")
        self.par_setd("domain1", "x2max", "%.15e", ylength, "-")
        self.par_setd("domain1", "x3min", "%.15e", 0.0, "-")
        self.par_setd("domain1", "x3max", "%.15e", zlength, "-")
        self.par_seti("domain1", "iDisp", "%d", 0, "-")
        self.par_seti("domain1", "jDisp", "%d", 0, "-")
        self.par_seti("domain1", "kDisp", "%d", 0, "-")
        
        return 0
    
    def define_subgrid(self, level, nmeshx, nmeshy, nmeshz, i, j, k):
        """
        Define a new domain on the given level the number of cells in this 
        domain is given by nmeshx,  nmeshy, nmeshz. 
        
        Each level is twice as dense as in every directory as 
        the previous level (there are 8 cells per higher level cell).
        """
        self.number_of_grids += 1
        
        domain = "domain{0}".format(self.number_of_grids)
        
        self.par_seti("job","num_domains", "%d", self.number_of_grids, "-")
        
        self.par_seti(domain, "level", "%d", level, "-")
        self.par_seti(domain, "Nx1", "%d", nmeshx, "-")
        self.par_seti(domain, "Nx2", "%d", nmeshy, "-")
        self.par_seti(domain, "Nx3", "%d", nmeshz, "-")
        self.par_seti(domain, "iDisp", "%d", i, "-")
        self.par_seti(domain, "jDisp", "%d", j, "-")
        self.par_seti(domain, "kDisp", "%d", k, "-")
        
        return self.number_of_grids

    @legacy_function    
    def get_position_of_index():
        """
        Retrieves the x, y and z position of the center of
        the cell with coordinates i, j, k in the grid specified
        by the index_of_grid
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'i'
        return function
    
    '''
    @legacy_function    
    def get_index_of_position():
        """
        Retrieves the i,j and k index of the grid cell containing the
        given x, y and z position. The cell is looked up
        in the grid specified by index_of_grid.
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        
        for x in ['i','j','k']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)

        function.result_type = 'i'
        return function
    '''
    def get_index_range_inclusive(self, index_of_grid = 1):
        """
        Returns the min and max values of indices in each
        direction. The range is inclusive, the min index
        and max index both exist and can be queried.
        The total number of cells in one direction
        is max - min + 1.
        """
        domainid = "domain{0}".format(index_of_grid)
        
        ni = self.par_geti(domainid, "Nx1")
        nj = self.par_geti(domainid, "Nx2")
        nk = self.par_geti(domainid, "Nx3")
        idisp = self.par_geti(domainid, "iDisp")
        jdisp = self.par_geti(domainid, "jDisp")
        kdisp = self.par_geti(domainid, "kDisp")
        #print index_of_grid, "  ===  > ", (idisp, idisp+ni[0]-1, jdisp, jdisp + nj[0]-1, kdisp, kdisp + nk[0]-1)
        
        return (idisp, idisp+ni-1, jdisp, jdisp + nj-1, kdisp, kdisp + nk-1)
        
    

    def get_index_range_magnetic_field_inclusive(self, index_of_grid = 1):
        original = list(self.get_index_range_inclusive(index_of_grid))
        original[1] += 1
        original[3] += 1
        original[5] += 1
        return original
        
    def get_mesh_indices(self):
        """
        Return 3 arrays, containing the indices for i, j and k
        """
        si,ei,sj,ej,sk,ek = self.get_index_range_inclusive()
        indexgrid = numpy.mgrid[slice(si,ei+1),slice(sj,ej+1),slice(sk,ek+1)]
        return indexgrid.reshape(3, -1)
        
    
        
    def set_four_pi_G(self, value):
        self.par_setd("problem", "four_pi_G", "%.15e", value, "")
        return 0 
        
    def set_grav_mean_rho(self, value):
        self.par_setd("problem", "grav_mean_rho", "%.15e", value, "")
        return 0 

    def set_isocsound(self, value):
        self.par_setd("problem", "iso_csound", "%.15e", value, "")
        return 0 
        
    def set_gamma(self, value):
        self.par_setd("problem", "gamma", "%.15e", value, "-") 
        return 0 
    
    def set_courant_friedrichs_lewy_number(self, value):
        self.par_setd("time", "cour_no", "%.15e", value, "-")
        return 0 
        
    def set_boundary(self, xbound1, xbound2, ybound1, ybound2, zbound1, zbound2):
        map_from_string_to_flag = {
            "reflective": 1, 
            "outflow":2, 
            "periodic":4,
            "interface": 10,
        }
        
        self.par_seti("domain1", "bc_ix1", "%d", map_from_string_to_flag[xbound1], "-")
        self.par_seti("domain1", "bc_ox1", "%d", map_from_string_to_flag[xbound2], "-")
        self.par_seti("domain1", "bc_ix2", "%d", map_from_string_to_flag[ybound1], "-")
        self.par_seti("domain1", "bc_ox2", "%d", map_from_string_to_flag[ybound2], "-")
        self.par_seti("domain1", "bc_ix3", "%d", map_from_string_to_flag[zbound1], "-")
        self.par_seti("domain1", "bc_ox3", "%d", map_from_string_to_flag[zbound2], "-")
        
        return 0
        
    
    def set_parallel_decomposition(self, nx, ny, nz):
        if nx == 0 or ny == 0 or nz == 0:
            self.par_seti("domain1", "AutoWithNProc", "%d", self.channel.number_of_workers, "-")
        else:
            self.par_seti("domain1", "AutoWithNProc", "%d", 0, "-")
            self.par_seti("domain1", "NGrid_x1", "%d", nx, "-")
            self.par_seti("domain1", "NGrid_x2", "%d", ny, "-")
            self.par_seti("domain1", "NGrid_x3", "%d", nz, "-")
        return 0
        
    def set_auto_decomposition(self, value):
        self.par_seti("parallel", "auto", "%d", value, "-")
        return 0
        
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
    def get_evolve_to_exact_time():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='bool', direction=function.OUT) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_evolve_to_exact_time():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='bool', direction=function.IN) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_has_external_gravitational_potential():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='int32', direction=function.IN) 
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_has_external_gravitational_potential():
        function = LegacyFunctionSpecification() 
        function.addParameter('value', dtype='int32', direction=function.OUT) 
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
        
    def get_index_range_for_potential(self, index_of_grid = 1):
        """
        Returns the min and max values of indices in each
        direction for the potential field, this
        range is 1 cell larger than the normal grid
        in all directions"""
        imin,imax,jmin,jmax,kmin,kmax = numpy.asarray(self.get_index_range_inclusive(index_of_grid = index_of_grid))
        imin -= 1
        imax += 1
        if jmin != jmax:
            jmin -= 1
            jmax += 1
        if kmin != kmax:
            kmin -= 1
            kmax += 1
        return imin,imax, jmin, jmax, kmin, kmax        
    
    @legacy_function    
    def set_potential():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        
        #function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        
        function.addParameter('potential', dtype='d', direction=function.IN)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_potential():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        #function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('potential', dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_interpolated_gravitational_potential():
        """
        Return the interpolated gravitational potential, can
        only interpolate over one axis at the time and
        only at half way points between the grid points.
        **For debugging purposes only**
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('potential', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    def get_isocsound(self):
        return self.par_getd("problem", "iso_csound"), 0

    def get_gamma(self):
        return self.par_getd("problem", "gamma"), 0
        
    def get_four_pi_G(self):
        return self.par_getd("problem", "four_pi_G"), 0
        
    def get_courant_friedrichs_lewy_number(self):
        return self.par_getd("time", "cour_no"), 0
        
    def get_grav_mean_rho(self):
        return self.par_getd("problem", "grav_mean_rho"), 0
        
    @legacy_function
    def get_grid_gravitational_potential():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['phi']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        
        return function
    
    
    

    @legacy_function
    def get_grid_acceleration():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['fx', 'fy', 'fz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        
        return function
        
    @legacy_function
    def set_grid_acceleration():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['fx', 'fy', 'fz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        
        return function

    @legacy_function
    def get_grid_gravitational_acceleration():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['fx', 'fy', 'fz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        
        return function
        
    @legacy_function
    def get_gravity_at_point():
        """
        Determine the gravitational force on a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eps', dtype='float64', direction=function.IN,
            description = "The smoothing parameter")
        function.addParameter('x', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('y', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('z', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('forcex', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcey', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcez', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
         0 - OK
            Force could be calculated
        -1 - ERROR
            No force calculation supported
        """
        return function


    @legacy_function
    def get_potential_at_point():
        """
        Determine the potential on a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eps', dtype='float64', direction=function.IN,
         description = "The smoothing factor, may be ignored by the code")
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('phi', dtype='float64', direction=function.OUT)
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
        

    @legacy_function
    def get_grid_scalar():
        """
        Retreives advected scalar property
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['scalar',]:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def set_grid_scalar():
        """
        Stores advected scalar property
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['scalar',]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    
    
    @legacy_function
    def get_boundary_index_range_inclusive():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
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
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
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
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['rho','rhovx','rhovy','rhovz','en']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    
    @legacy_function    
    def get_boundary_position_of_index():
        """
        Retrieves the x, y and z position of the center of
        the cell with coordinates i, j, k in the grid specified
        by the index_of_grid
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN, default = 1)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'i'
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
        
    @legacy_function    
    def get_hydro_state_for_cell():
        function = LegacyFunctionSpecification()  
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['dx','dy','dz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, default = 0)
        for x in ['rho','rhovx','rhovy','rhovz','rhoe']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'i' 
        function.must_handle_array = True
        return function
    
class Athena(CommonCode):

    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        
        self.stopping_conditions = StoppingConditions(self)
        CommonCode.__init__(self,  AthenaInterface(**options), **options)
        
    def define_converter(self, handler):
        if self.unit_converter is None:
            return
        
        handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())

    def define_properties(self, handler):
        handler.add_property('get_time', public_name = "model_time")
        
    def define_methods(self, handler):
        handler.add_method(
            'evolve_model',
            (time,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_position_of_index',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (length, length, length, handler.ERROR_CODE,)
        )
        
        density = mass / (length**3)
        momentum =  mass / (time * (length**2))
        energy =  mass / ((time**2) * length)
        potential_energy =  length ** 2 / time ** 2
        magnetic_field = mass / current / time ** 2
        
        handler.add_method(
            'set_grid_state',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            density, momentum, momentum, momentum, energy,
            handler.INDEX),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'set_grid_magnetic_field',
            (handler.INDEX, handler.INDEX, handler.INDEX,
             magnetic_field, magnetic_field, magnetic_field,
            handler.INDEX),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            'get_grid_state',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (density, momentum, momentum, momentum, energy,
            handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grid_energy_density',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            energy, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_grid_energy_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            ( energy,
            handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grid_density',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            density, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_grid_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (density,
            handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grid_scalar',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            handler.NO_UNIT, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_grid_scalar',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (handler.NO_UNIT,
            handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grid_momentum_density',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            momentum, momentum, momentum, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_grid_momentum_density',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            ( momentum, momentum, momentum, 
            handler.ERROR_CODE,)
        )
    
        handler.add_method(
            'get_grid_magnetic_field',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            ( magnetic_field, magnetic_field, magnetic_field,
            handler.ERROR_CODE,)
        )

        handler.add_method(
            'set_potential',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            potential_energy),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_potential',
            (handler.INDEX, handler.INDEX, handler.INDEX,),
            (potential_energy,
            handler.ERROR_CODE,)
        )
    
        handler.add_method(
            'get_grid_gravitational_potential',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX,),
            (potential_energy,
            handler.ERROR_CODE,)
        )
         
        handler.add_method(
            'get_interpolated_gravitational_potential',
            (length, length, length),
            (potential_energy,
            handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_grid_gravitational_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX,),
            (acceleration,acceleration,acceleration,
            handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_grid_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX,),
            (acceleration,acceleration,acceleration,
            handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grid_acceleration',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            acceleration,acceleration,acceleration, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_gravity_at_point',
            (length, length, length, length),
            (acceleration, acceleration, acceleration, handler.ERROR_CODE)
        )

        handler.add_method(
            'get_potential_at_point',
            (length, length, length, length),
            (potential, handler.ERROR_CODE)
        )
    
        handler.add_method(
            "get_isocsound",
            (),
            (length / time, handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "set_isocsound",
            (length / time, ),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_timestep",
            (),
            (time, handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "set_timestep",
            (time, ),
            (handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "get_gamma",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "set_gamma",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "get_courant_friedrichs_lewy_number",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
    
        handler.add_method(
            "set_courant_friedrichs_lewy_number",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
    
        handler.add_method(
            'get_time',
            (),
            (time, handler.ERROR_CODE,)
        )
    
        handler.add_method(
            'set_four_pi_G',
            ( length**3 / (mass * time**2)),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_four_pi_G',
            (),
            ( (length**3) / (mass * (time**2)), handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_grav_mean_rho',
            (  mass / length**3, ),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_grav_mean_rho',
            (),
            (mass / length**3, handler.ERROR_CODE,)
        )
    
        handler.add_method(
            'setup_mesh',
            (handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, length, length, length,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            'set_boundary',
            (handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'set_boundary_state',
            (handler.INDEX, handler.INDEX, handler.INDEX,
            density, momentum, momentum, momentum, energy,
            handler.INDEX, handler.INDEX),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            'get_boundary_state',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (density, momentum, momentum, momentum, energy,
            handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_boundary_position_of_index',
            (handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX, handler.INDEX),
            (length, length, length, handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_boundary_index_range_inclusive',
            (handler.INDEX, handler.INDEX),
            (handler.NO_UNIT, handler.NO_UNIT,handler.NO_UNIT, handler.NO_UNIT,handler.NO_UNIT, handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            'get_hydro_state_at_point',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,
                generic_unit_system.speed, generic_unit_system.speed, generic_unit_system.speed),
            (generic_unit_system.density, generic_unit_system.momentum_density, generic_unit_system.momentum_density, 
                generic_unit_system.momentum_density, generic_unit_system.energy_density, handler.ERROR_CODE)
        )
        handler.add_method(
            'get_hydro_state_for_cell',
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,
                generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,
                generic_unit_system.speed, generic_unit_system.speed, generic_unit_system.speed),
            (generic_unit_system.density, generic_unit_system.momentum_density, generic_unit_system.momentum_density, 
                generic_unit_system.momentum_density, generic_unit_system.energy_density, handler.ERROR_CODE)
        )
        
        
        
        self.stopping_conditions.define_methods(handler)
    
    
    def specify_grid(self, definition, index_of_grid = 1):
        definition.set_grid_range('get_index_range_inclusive')
        
        definition.add_getter('get_position_of_index', names=('x','y','z'))
        
        definition.add_getter('get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        
        definition.add_getter('get_grid_density', names=('rho',))
        definition.add_setter('set_grid_density', names=('rho',))
        
        if self.mode == self.MODE_SCALAR:
            definition.add_getter('get_grid_scalar', names=('scalar',))
            definition.add_setter('set_grid_scalar', names=('scalar',))
            
        definition.add_getter('get_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_setter('set_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_getter('get_grid_energy_density', names=('energy',))
        definition.add_setter('set_grid_energy_density', names=('energy',))
        
        
        definition.add_getter('get_grid_gravitational_potential', names=('gravitational_potential',))
        definition.add_getter('get_grid_gravitational_acceleration', names=('gravitational_acceleration_x','gravitational_acceleration_y','gravitational_acceleration_z',))
        
        definition.add_getter('get_grid_acceleration', names=('ax','ay','az'))
        definition.add_setter('set_grid_acceleration', names=('ax','ay','az'))
        
        definition.define_extra_keywords({'index_of_grid':index_of_grid})
        

    def specify_mangnetic_filed_grid(self, definition, index_of_grid = 1):
        definition.set_grid_range('get_index_range_magnetic_field_inclusive')
        
        definition.add_getter('get_position_of_index', names=('x','y','z'))
    
        definition.add_getter('get_grid_magnetic_field', names=('B1i','B2i','B3i'))   
        definition.add_setter('set_grid_magnetic_field', names=('B1i','B2i','B3i'))
         
        definition.define_extra_keywords({'index_of_grid':index_of_grid})
        
  
    def specify_boundary_grid(self, definition, index_of_boundary, index_of_grid = 1):
        definition.set_grid_range('get_boundary_index_range_inclusive')
        
        definition.add_getter('get_boundary_position_of_index', names=('x','y','z'))
        
        definition.add_getter('get_boundary_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_boundary_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
       
        definition.define_extra_keywords({'index_of_boundary': index_of_boundary, 'index_of_grid':index_of_grid})
        
    
    
    def sepecify_extended_grid(self, definition, index_of_grid = 1):
        self.specify_grid(definition, index_of_grid = index_of_grid)
        definition.set_grid_range('get_index_range_extended')
        
        
    @property
    def grid(self):
        return self._create_new_grid(self.specify_grid, index_of_grid = 1)
    
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
    

    def get_extended_grid(self, index_of_grid = 1):
        return self._create_new_grid(self.sepecify_extended_grid, index_of_grid = index_of_grid)
                
    def get_index_range_extended(self, index_of_grid = 1):
        i0,i1, j0,j1, k0,k1 = self.get_index_range_inclusive(index_of_grid = index_of_grid)
        dj = 2 if j1 > j0 else 0
        dk = 2 if k1 > k0 else 0
        return i0-2, i1+2, j0-dj, j0+dj, k0-dk, k1+dk
    
    def itergrids(self):
        n = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield self._create_new_grid(self.specify_grid, index_of_grid = x)

    def iter_magnetic_field_grids(self):
        n = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield self._create_new_grid(self.specify_mangnetic_filed_grid, index_of_grid = x)

    
    def iter_hydro_and_mhd_grids(self):
        n = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield (
                self._create_new_grid(self.specify_grid, index_of_grid = x),
                self._create_new_grid(self.specify_mangnetic_filed_grid, index_of_grid = x),
            )

    def define_particle_sets(self, handler):
        handler.define_grid('potential_grid')
        handler.set_grid_range('potential_grid', 'get_index_range_for_potential')
        handler.add_getter('potential_grid', 'get_position_of_index', names=('x','y','z'))
        handler.add_getter('potential_grid', 'get_potential', names=('potential',))
        handler.add_setter('potential_grid', 'set_potential', names=('potential', ))
        handler.define_extra_keywords('potential_grid', {'index_of_grid':1})
        
        
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_isocsound", 
            "set_isocsound",
            "isothermal_sound_speed", 
            "isothermal sound speed, only used for isothermal EOS", 
            default_value = 0.0 | length / time,
            must_set_before_get = True
        )
        
        handler.add_method_parameter(
            "get_gamma", 
            "set_gamma",
            "gamma", 
            "ratio of specific heats used in equation of state", 
            default_value = 1.6666666666666667,
            must_set_before_get = True
        )
        
        handler.add_method_parameter(
            "get_four_pi_G", 
            "set_four_pi_G",
            "four_pi_G", 
            "value of four times pi time G", 
            default_value = 4 * numpy.pi * (1| (length**3) / (mass * (time**2))),
            must_set_before_get = True
        )
        
        handler.add_method_parameter(
            "get_grav_mean_rho", 
            "set_grav_mean_rho",
            "gravity_mean_rho", 
            "define the mean density in the field for self gravity calulations", 
            default_value = 0 | mass / length ** 3,
            must_set_before_get = True
        )
        
        
        handler.add_method_parameter(
            "get_courant_friedrichs_lewy_number", 
            "set_courant_friedrichs_lewy_number",
            "courant_number", 
            "CFL number", 
            default_value = 0.3,
            must_set_before_get = True
        )
        
        
        handler.add_method_parameter(
            "get_evolve_to_exact_time", 
            "set_evolve_to_exact_time",
            "must_evolve_to_exact_time", 
            "End the evolve model at the exact specified time", 
            default_value = True
        )
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "nmeshx",
            "nx", 
            "number of cells in the x direction", 
            10,
        )
        
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "nmeshy",
            "ny", 
            "number of cells in the y direction", 
            10,
        )
        
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "nmeshz",
            "nz", 
            "number of cells in the z direction", 
            10,
        )
        
        handler.add_caching_parameter(
            "setup_mesh", 
            "xlength",
            "length_x", 
            "length of model in the x direction", 
            10 | length,
        )
        handler.add_caching_parameter(
            "setup_mesh", 
            "ylength",
            "length_y", 
            "length of model in the x direction", 
            10 | length,
        )
        handler.add_caching_parameter(
            "setup_mesh", 
            "zlength",
            "length_z", 
            "length of model in the z direction", 
            10 | length,
        )
        
        handler.add_vector_parameter(
            "mesh_size",
            "number of cells in the x, y and z directions",
            ("nx", "ny", "nz")
        )
        
        handler.add_vector_parameter(
            "mesh_length",
            "length of the model in the x, y and z directions",
            ("length_x", "length_y", "length_z")
        )
        
        handler.add_caching_parameter(
            "set_parallel_decomposition", 
            "nx",
            "nproc_x", 
            "number of processors for the x direction",
            0,
        )
        handler.add_caching_parameter(
            "set_parallel_decomposition", 
            "ny",
            "nproc_y", 
            "number of processors for the y direction",
            0,
        )
        handler.add_caching_parameter(
            "set_parallel_decomposition", 
            "nz",
            "nproc_z", 
            "number of processors for the z direction",
            0,
        )
        
        handler.add_vector_parameter(
            "parallel_decomposition",
            "number of processors for each dimensions",
            ("nproc_x", "nproc_y", "nproc_z")
        )
        
        handler.add_caching_parameter(
            "set_boundary", 
            "xbound1",
            "xbound1", 
            "boundary conditions on first (inner, left) X boundary", 
            "reflective",
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "xbound2",
            "xbound2", 
            "boundary conditions on second (outer, right) X boundary",
            "reflective",
        )
        
        handler.add_caching_parameter(
            "set_boundary", 
            "ybound1",
            "ybound1", 
            "boundary conditions on first (inner, front) Y boundary", 
            "reflective",
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "ybound2",
            "ybound2", 
            "boundary conditions on second (outer, back) Y boundary",
            "reflective",
        )
        
        handler.add_caching_parameter(
            "set_boundary", 
            "zbound1",
            "zbound1", 
            "boundary conditions on first (inner, bottom) Z boundary", 
            "reflective",
        )
        
        
        handler.add_caching_parameter(
            "set_boundary", 
            "zbound2",
            "zbound2", 
            "boundary conditions on second (outer, top) Z boundary", 
            "reflective",
        )
        
        
        
        handler.add_vector_parameter(
            "x_boundary_conditions",
            "boundary conditions for the X directorion",
            ("xbound1", "xbound2")
        )
        
        
        handler.add_vector_parameter(
            "y_boundary_conditions",
            "boundary conditions for the Y directorion",
            ("ybound1", "ybound2")
        )
        
        
        handler.add_vector_parameter(
            "z_boundary_conditions",
            "boundary conditions for the Z directorion",
            ("zbound1", "zbound2")
        )
        
        self.stopping_conditions.define_parameters(handler)

    def commit_parameters(self):
        self.parameters.send_not_set_parameters_to_code()
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
    
    def define_state(self, handler): 
        CommonCode.define_state(self, handler)       
        #handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        
        handler.add_transition('EDIT', 'RUN', 'initialize_grid')
        handler.add_method('RUN', 'evolve_model')
        handler.add_method('RUN', 'get_hydro_state_at_point')
        
        for state in ['EDIT', 'RUN']:
            for methodname in [
                    'get_grid_state',
                    'set_grid_state',
                    'get_potential',
                    'set_potential',
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
                    'get_number_of_grids',
                    'get_index_range_inclusive',
                    'get_boundary_state',
                    'set_boundary_state',
                    'get_boundary_position_if_index',
                    'get_boundary_index_range_inclusive'
                ]:
                handler.add_method(state, methodname)
                
        self.stopping_conditions.define_state(handler)

