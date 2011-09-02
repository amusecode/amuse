from amuse.community import *



from amuse.community.interface.mhd import MagnetohydrodynamicsInterface

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
        .. [#] Gardiner & Stone 2005, JCP, 205, 509  (2D JCP Method)
        .. [#] Gardiner & Stone 2007, JCP, 227, 4123 (3D JCP Method)
        .. [#] Stone et al. 2008, ApJS, 178, 137 (Method)
        .. [#] Stone & Gardiner 2009, NewA, 14, 139 (van Leer Integrator)
        .. [#] Skinner & Ostriker 2010, ApJ, 188, 290 (Cylindrical Integrator)
        .. [#] Stone & Gardiner 2010, ApJS, 189, 142 (Shearing Box Method)
    """
    
    include_headers = ['worker_code.h', 'stopcond.h']
    
    MODE_NORMAL = 'normal'
    MODE_SELF_GRAVITY   = 'self-gravity'
    MODE_MHD   = 'mhd'
    
    def __init__(self, mode = MODE_NORMAL, **options):
        
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **options)
        self.set_auto_decomposition(1)
        LiteratureReferencesMixIn.__init__(self)
        self.number_of_grids = 1
        
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'athena_worker'
        elif mode == self.MODE_SELF_GRAVITY:
            return 'athena_worker_selfgrav'
        elif mode == self.MODE_MHD:
            return 'athena_worker_mhd'
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
        self.par_seti("domain1", "AutoWithNProc", "%d", self.channel.number_of_workers, "-")
        
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
        idisp = self.par_geti(domainid, "iDisp")[0]
        jdisp = self.par_geti(domainid, "jDisp")[0]
        kdisp = self.par_geti(domainid, "kDisp")[0]
        #print index_of_grid, "  ===  > ", (idisp, idisp+ni[0]-1, jdisp, jdisp + nj[0]-1, kdisp, kdisp + nk[0]-1)
        
        return (idisp, idisp+ni[0]-1, jdisp, jdisp + nj[0]-1, kdisp, kdisp + nk[0]-1)

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
        map_from_string_to_flag = {"reflective": 1, "outflow":2, "periodic":4}
        
        self.par_seti("domain1", "bc_ix1", "%d", map_from_string_to_flag[xbound1], "-")
        self.par_seti("domain1", "bc_ox1", "%d", map_from_string_to_flag[xbound2], "-")
        self.par_seti("domain1", "bc_ix2", "%d", map_from_string_to_flag[ybound1], "-")
        self.par_seti("domain1", "bc_ox2", "%d", map_from_string_to_flag[ybound2], "-")
        self.par_seti("domain1", "bc_ix3", "%d", map_from_string_to_flag[zbound1], "-")
        self.par_seti("domain1", "bc_ox3", "%d", map_from_string_to_flag[zbound1], "-")
        
        return 0
        
    
    def set_parallel(self, nx, ny, nz):
        self.par_seti("parallel", "NGrid_x1", "%d", nx, "-")
        self.par_seti("parallel", "NGrid_x2", "%d", ny, "-")
        self.par_seti("parallel", "NGrid_x3", "%d", nz, "-")
        
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
            kmin -=1
            kmax +=1
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
    
    def get_courant_friedrichs_lewy_number(self):
        return self.par_getd("time", "cour_no"), 0
    

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
        

    
    
    
class Athena(InCodeComponentImplementation):

    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        
        self.stopping_conditions = StoppingConditions(self)
        InCodeComponentImplementation.__init__(self,  AthenaInterface(**options), **options)
        
    def define_converter(self, object):
        if self.unit_converter is None:
            return
        
        object.set_converter(self.unit_converter.as_converter_from_si_to_generic())

    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")
        
    def define_methods(self, object):
        object.add_method(
            'evolve_model',
            (time,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length, length, length, object.ERROR_CODE,)
        )
        
        density = mass / (length**3)
        momentum =  mass / (time * (length**2))
        energy =  mass / ((time**2) * length)
        potential_energy =  length ** 2 / time ** 2
        magnetic_field = mass / current / time ** 2
        
        object.add_method(
            'set_grid_state',
            (object.INDEX, object.INDEX, object.INDEX,
            density, momentum, momentum, momentum, energy,
            object.INDEX),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            'set_grid_magnetic_field',
            (object.INDEX, object.INDEX, object.INDEX,
             magnetic_field, magnetic_field, magnetic_field,
            object.INDEX),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density, momentum, momentum, momentum, energy,
            object.ERROR_CODE,)
        )
        object.add_method(
            'get_grid_energy_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            ( energy,
            object.ERROR_CODE,)
        )
        object.add_method(
            'get_grid_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density,
            object.ERROR_CODE,)
        )
        object.add_method(
            'get_grid_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            ( momentum, momentum, momentum, 
            object.ERROR_CODE,)
        )
    
        object.add_method(
            'get_grid_magnetic_field',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            ( magnetic_field, magnetic_field, magnetic_field,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_potential',
            (object.INDEX, object.INDEX, object.INDEX,
            potential_energy),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'get_potential',
            (object.INDEX, object.INDEX, object.INDEX,),
            (potential_energy,
            object.ERROR_CODE,)
        )
    
        object.add_method(
            'get_grid_gravitational_potential',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX,),
            (potential_energy,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_grid_gravitational_acceleration',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX,),
            (acceleration,acceleration,acceleration,
            object.ERROR_CODE,)
        )
        
        object.add_method(
            'get_gravity_at_point',
            (length, length, length, length),
            (acceleration, acceleration, acceleration, object.ERROR_CODE)
        )

        object.add_method(
            'get_potential_at_point',
            (length, length, length, length),
            (potential, object.ERROR_CODE)
        )
    
        object.add_method(
            "get_isocsound",
            (),
            (length / time, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_isocsound",
            (length / time, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            "get_gamma",
            (),
            (units.none, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_gamma",
            (units.none, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            "get_courant_friedrichs_lewy_number",
            (),
            (units.none, object.ERROR_CODE,)
        )
    
        object.add_method(
            "set_courant_friedrichs_lewy_number",
            (units.none, ),
            (object.ERROR_CODE,)
        )
    
        object.add_method(
            'get_time',
            (),
            (time, object.ERROR_CODE,)
        )
    
        object.add_method(
            'setup_mesh',
            (units.none, units.none, units.none, length, length, length,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'set_boundary',
            (units.string, units.string, units.string, units.string, units.string, units.string,),
            (object.ERROR_CODE,)
        )
        self.stopping_conditions.define_methods(object)
    
    
    def specify_grid(self, definition, index_of_grid = 1):
        definition.set_grid_range('get_index_range_inclusive')
        
        definition.add_getter('get_position_of_index', names=('x','y','z'))
        
        definition.add_getter('get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        
        definition.add_getter('get_grid_density', names=('rho',))
        definition.add_getter('get_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_getter('get_grid_energy_density', names=('energy',))
        
        
        definition.add_getter('get_grid_gravitational_potential', names=('gravitational_potential',))
        definition.add_getter('get_grid_gravitational_acceleration', names=('gravitational_acceleration_x','gravitational_acceleration_y','gravitational_acceleration_z',))
        
        definition.define_extra_keywords({'index_of_grid':index_of_grid})
        

    def specify_mangnetic_filed_grid(self, definition, index_of_grid = 1):
        definition.set_grid_range('get_index_range_magnetic_field_inclusive')
        
        definition.add_getter('get_position_of_index', names=('x','y','z'))
    
        definition.add_getter('get_grid_magnetic_field', names=('B1i','B2i','B3i'))   
        definition.add_setter('set_grid_magnetic_field', names=('B1i','B2i','B3i'))
         
        definition.define_extra_keywords({'index_of_grid':index_of_grid})
        

    @property
    def grid(self):
        return self._create_new_grid(self.specify_grid, index_of_grid = 1)
        
    def itergrids(self):
        n, error = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield self._create_new_grid(self.specify_grid, index_of_grid = x)

    def iter_magnetic_field_grids(self):
        n, error = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield self._create_new_grid(self.specify_mangnetic_filed_grid, index_of_grid = x)

    
    def iter_hydro_and_mhd_grids(self):
        n, error = self.get_number_of_grids()
        
        for x in range(1,n+1):
            yield (
                self._create_new_grid(self.specify_grid, index_of_grid = x),
                self._create_new_grid(self.specify_mangnetic_filed_grid, index_of_grid = x),
            )

    def define_particle_sets(self, object):
        object.define_grid('potential_grid')
        object.set_grid_range('potential_grid', 'get_index_range_for_potential')
        object.add_getter('potential_grid', 'get_position_of_index', names=('x','y','z'))
        object.add_getter('potential_grid', 'get_potential', names=('potential',))
        object.add_setter('potential_grid', 'set_potential', names=('potential', ))
        object.define_extra_keywords('potential_grid', {'index_of_grid':1})
        
        
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_isocsound", 
            "set_isocsound",
            "isothermal_sound_speed", 
            "isothermal sound speed, only used for isothermal EOS", 
            default_value = 0.0 | length / time,
            must_set_before_get = True
        )
        
        object.add_method_parameter(
            "get_gamma", 
            "set_gamma",
            "gamma", 
            "ratio of specific heats used in equation of state", 
            default_value = 1.6666666666666667 | units.none,
            must_set_before_get = True
        )
        
        object.add_method_parameter(
            "get_courant_friedrichs_lewy_number", 
            "set_courant_friedrichs_lewy_number",
            "courant_number", 
            "CFL number", 
            default_value = 0.3 | units.none,
            must_set_before_get = True
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshx",
            "nx", 
            "number of cells in the x direction", 
            10 | units.none,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshy",
            "ny", 
            "number of cells in the y direction", 
            10 | units.none,
        )
        
        
        object.add_caching_parameter(
            "setup_mesh", 
            "nmeshz",
            "nz", 
            "number of cells in the z direction", 
            10 | units.none,
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
            "set_boundary", 
            "xbound1",
            "xbound1", 
            "boundary conditions on first (inner, left) X boundary", 
            "reflective" | units.string,
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "xbound2",
            "xbound2", 
            "boundary conditions on second (outer, right) X boundary",
            "reflective" | units.string,
        )
        
        object.add_caching_parameter(
            "set_boundary", 
            "ybound1",
            "ybound1", 
            "boundary conditions on first (inner, front) Y boundary", 
            "reflective" | units.string,
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "ybound2",
            "ybound2", 
            "boundary conditions on second (outer, back) Y boundary",
            "reflective" | units.string,
        )
        
        object.add_caching_parameter(
            "set_boundary", 
            "zbound1",
            "zbound1", 
            "boundary conditions on first (inner, bottom) Z boundary", 
            "reflective" | units.string,
        )
        
        
        object.add_caching_parameter(
            "set_boundary", 
            "zbound2",
            "zbound2", 
            "boundary conditions on second (outer, top) Z boundary", 
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
        
        self.stopping_conditions.define_parameters(object)

    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
    
    
