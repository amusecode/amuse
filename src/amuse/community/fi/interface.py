import os
import numpy
from amuse.community.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
from amuse.community import *
from amuse.support.options import option

class FiInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn, StoppingConditionInterface):   
    """
    FI is a parallel TreeSPH code for galaxy simulations. Extensively 
    rewritten, extended and parallelized it is a development from code from 
    Jeroen Gerritsen and Roelof Bottema, which itself goes back to Treesph. 
    
    Note that some features are not working atm. These may be fixed in the 
    future. (and I will think of a better name)
    
    The relevant references are:
        .. [#] Pelupessy, PhD thesis 2005, Leiden Observatory
        .. [#] Pelupessy, van der Werf & Icke, 2004, A&A 422, 55
        .. [#] Gerritsen & Icke, 1997, A&A 325, 972
        .. [#] Hernquist & Katz, 1989, ApJS 70, 419
    """
    get_total_radius=None
    get_total_mass=None
    get_center_of_mass_position=None
    get_center_of_mass_velocity=None
    set_acceleration=None
    get_acceleration=None
    
    
    MODE_NORMAL = 'normal'
    MODE_PERIODIC_BOUNDARIES   = 'periodic'
    
    def __init__(self, mode = MODE_NORMAL,  **options):
        self.mode = mode
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
                     
    
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'fi_worker'
        elif mode == self.MODE_PERIODIC_BOUNDARIES:
            return 'fi_worker_periodic'
        else:
            return 'fi_worker'
        
            
    @option(type="string", sections=('data',))
    def input_data_root_directory(self):
        """
        The root directory of the input data, read only directories
        """
        return os.path.join(get_amuse_root_dir(), 'data')
        
    @option(type="string", sections=('data',))
    def output_data_root_directory(self):
        """
        The root directory of the output data,
        read - write directory
        """
        return os.path.join(get_amuse_root_dir(), 'data')
        
    def get_data_directory(self):
        """
        Returns the root name of the directory for the Fi
        application data files.
        """
        return os.path.join(self.input_data_root_directory, 'fi', 'input')
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, 'data', 'fi', 'output')
    
    
        
    
       
    def new_particle(self, mass, x, y, z, vx, vy, vz, radius = 0.0):
        return self.new_dm_particle(mass, x, y, z, vx, vy, vz, radius)
    
    @legacy_function    
    def new_dm_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.IN, default = 0)
            
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_sph_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('h_smooth', dtype='d', direction=function.IN, default = 0)
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_star_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('tform', dtype='d', direction=function.IN, default = 0)
        function.addParameter('radius', dtype='d', direction=function.IN, default = 0)
        function.result_type = 'i'
        return function
             
    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz','radius']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state_sph():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz','u','h_smooth']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state_star():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz','tform','radius']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.IN, default = 0)
            
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state_star():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        
        function.addParameter('tform', dtype='d', direction=function.IN, default = 0)
        function.addParameter('radius', dtype='d', direction=function.IN, default = 0)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def set_state_sph():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','x','y','z','vx','vy','vz','u',]:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('h_smooth', dtype='d', direction=function.IN, default = 0)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_radius():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_position():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_velocity():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_internal_energy():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('u', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_dinternal_energy_dt():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('du_dt', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_smoothing_length():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('h_smooth', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('rho', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function    
    def get_star_tform():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('tform', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def set_mass():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_radius():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_position():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_velocity():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_internal_energy():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('u', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_smoothing_length():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('h_smooth', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def set_star_tform():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('tform', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    
    @legacy_function      
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('kinetic_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_potential_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_thermal_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('thermal_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_total_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('total_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function          
    def delete_particle():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function    

    @legacy_function    
    def evolve_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def get_number_of_particles():
        function = LegacyFunctionSpecification()  
        function.addParameter('number_of_particles', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function    
    def get_gravity_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['ax','ay','az']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'i' 
        function.must_handle_array = True
        return function

    @legacy_function    
    def get_potential_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('phi', dtype='d', direction=function.OUT)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'i' 
        function.must_handle_array = True
        return function

    @legacy_function    
    def get_hydro_state_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['rho','rhovx','rhovy','rhovz','rhoe']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'i' 
        function.must_handle_array = True
        return function

    @legacy_function    
    def synchronize_model():
        """ synchronize the model """
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

# setting/ getting parameters
# logicals
    @legacy_function   
    def set_use_hydro():
        """ set_use_hydro([0,1]): SPH hydro if 1, gravity only if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('use_hydro_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_use_hydro():
        """ get_use_hydro(): SPH hydro if 1, gravity only if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('use_hydro_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_radiate():
        """ set_radiate([0,1]): rad cooling/heating if 1, not if 0
             radiate false implies starform false """    
        function = LegacyFunctionSpecification()  
        function.addParameter('radiation_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_radiate():
        """ get_radiate(): rad cooling/heating if 1, not if 0
             radiate false implies starform false """        
        function = LegacyFunctionSpecification()  
        function.addParameter('radiation_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_starform():
        """ set_starform([0,1]): star formation if 1, not if 0 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('star_formation_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_starform():
        """ get_starform(): star formation if 1, not if 0 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('star_formation_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_cosmo():
        """ set_cosmo([0,1]): not functional at the moment """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_cosmo():
        """ get_cosmo(): not functional at the moment """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_sqrttstp():
        """ set_sqrttstp([0,1]): use sqrt(eps/acc) timestep crit if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('square_root_timestep_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sqrttstp():
        """ get_sqrttstp(): use sqrt(eps/acc) timestep crit if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('square_root_timestep_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_acc_tstp():
        """ set_acc_tstp([0,1]): use vref/acc timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('acc_timestep_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_acc_tstp():
        """ get_acc_tstp(): use vref/acc timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('acc_timestep_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_freetstp():
        """ set_freetstp([0,1]): use freeform (v/freev)**freevexp * (a/freea)**freeaexp timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('freeform_timestep_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_freetstp():
        """ get_freetstp(): use freeform timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('freeform_timestep_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_usequad():
        """ set_usequad([0,1]): calc. and use quadrupole cell moments if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('quadrupole_moments_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_usequad():
        """ get_usequad(): calc. and use quadrupole cell moments if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('quadrupole_moments_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_directsum():
        """ set_directsum([0,1]): direct N**2 grav sum if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('direct_sum_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_directsum():
        """ get_directsum(): direct N**2 grav sum if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('direct_sum_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_selfgrav():
        """ set_selfgrav([0,1]): calculate self-gravity if 1
          if set to 0, self gravity is not used, only external potentials"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('self_gravity_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_selfgrav():
        """ get_selfgrav(): calculate self-gravity if 1 """
        function = LegacyFunctionSpecification()  
        function.addParameter('self_gravity_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_fixthalo():
        """ set_fixthalo([0,1]): use fixed (spherical) potential if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('fixed_halo_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_fixthalo():
        """ get_fixthalo(): use fixed (spherical) potential if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('fixed_halo_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_adaptive_eps():
        """ set_adaptive_eps([0,1]): use of adaptive grav smoothing for all part if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('adaptive_smoothing_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_adaptive_eps():
        """ get_adaptive_eps(): use of adaptive grav smoothing for all part if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('adaptive_smoothing_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_gdgop():
        """ set_gdgop([0,1]): use of gadget cell opening criterion if 1 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_gdgop():
        """ get_gdgop(): use of gadget cell opening criterion if 1 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_smoothinput():
        """ set_smoothinput([0,1]): smooth input SPH prop. if 1 
         (not working) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('smooth_input_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_smoothinput():
        """ get_smoothinput(): smooth input SPH prop. if 1 
         (not working) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('smooth_input_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_consph():
        """ set_consph([0,1]): use springel&Hernquist conservative SPH form. if 1 
          at the moment this is only option"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('conservative_sph_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_consph():
        """ get_consph(): use springel&Hernquist conservative SPH form. if 1 
          at the moment this is only option"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('conservative_sph_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_sphinit():
        """ set_sphinit([0,1]): initialize sph dens and hsmooth if 1 
         most probably useless for AMUSE interface"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_dens_init_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sphinit():
        """ set_sphinit([0,1]): initialize sph dens and hsmooth if 1 
         most probably useless for AMUSE interface"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_dens_init_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_uentropy():
        """ set_uentropy([0,1]): integrate entropy if 1, internal energy if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('integrate_entropy_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_uentropy():
        """ get_uentropy(): integrate entropy if 1, internal energy if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('integrate_entropy_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_isotherm():
        """ set_isotherm([0,1]): isothermal gas if 1
          note that isotherm needs set_uentropy(0) (false)"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_isotherm():
        """ get_isotherm(): isothermal gas if 1
          note that isotherm needs set_uentropy(0) (false)"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_eps_is_h():
        """ set_eps_is_h([0,1]): gas particles grav. eps to SPH h if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('eps_is_h_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_eps_is_h():
        """ get_eps_is_h(): gas particles grav. eps to SPH h if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('eps_is_h_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;


# integers
    @legacy_function
    def set_firstsnap():
        """ no. of first snapshot """
        function = LegacyFunctionSpecification()  
        function.addParameter('first_snapshot', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_firstsnap():
        """ no. of first snapshot """
        function = LegacyFunctionSpecification()  
        function.addParameter('first_snapshot', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_stepout():
        """ no. of steps between output """
        function = LegacyFunctionSpecification()  
        function.addParameter('output_interval', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_stepout():
        """ no. of steps between output """
        function = LegacyFunctionSpecification()  
        function.addParameter('output_interval', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_steplog():
        """ no. of steps between logs """
        function = LegacyFunctionSpecification()  
        function.addParameter('log_interval', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_steplog():
        """ no. of steps between logs """
        function = LegacyFunctionSpecification()  
        function.addParameter('log_interval', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_max_tbin():
        """ maximum time bin (dtime*2**-max_tbin=minimum time step)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('maximum_time_bin', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_max_tbin():
        """ maximum time bin (dtime*2**-max_tbin=minimum time step)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('maximum_time_bin', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_minppbin():
        """ target no. of particles per time bin"""
        function = LegacyFunctionSpecification()  
        function.addParameter('minimum_part_per_bin', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_minppbin():
        """ target no. of particles per time bin"""
        function = LegacyFunctionSpecification()  
        function.addParameter('minimum_part_per_bin', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_targetnn():
        """ target no. of neighbour particles for variable grav. eps"""
        function = LegacyFunctionSpecification()  
        function.addParameter('targetnn', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_targetnn():
        """ target no. of neighbour particles for variable grav. eps"""
        function = LegacyFunctionSpecification()  
        function.addParameter('targetnn', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_verbosity():
        """ level of terminal output (0=minimum)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('verbosity', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_verbosity():
        """ level of terminal output (0=minimum)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('verbosity', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmooth', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmooth', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

# real
    @legacy_function
    def set_pboxsize():
        """ size of simulation domain box (particles outside get deleted)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('pboxsize', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_pboxsize():
        """ size of simulation domain box (particles outside get deleted)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('pboxsize', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_dtime():
        """ timestep (code units)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('dtime', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_dtime():
        """ timestep (code units)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('dtime', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_unitm_in_msun():
        """ code mass unit (in Msun, 10^9 standard) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitm_in_msun', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_unitm_in_msun():
        """ code mass unit (in Msun, 10^9 standard) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitm_in_msun', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_unitl_in_kpc():
        """ code length unit (in kpc, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitl_in_kpc', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_unitl_in_kpc():
        """ code length unit (in kpc, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitl_in_kpc', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;


    @legacy_function
    def set_tstepcrit():
        """ sqrttstp timestep constant (unitless,standard=1.) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstepcrit', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_tstepcrit():
        """ sqrttstp timestep constant (unitless,standard=1.) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstepcrit', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_tstpcr2():
        """ acc_tstp timestep constant (unitless,standard=0.25) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstpcr2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_tstpcr2():
        """ acc_tstp timestep constant (unitless,standard=0.25) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstpcr2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_freev():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freev', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_freev():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freev', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_freea():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freea', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_freea():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freea', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_freevexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freevexp', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_freevexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freevexp', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_freeaexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freeaexp', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_freeaexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freeaexp', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('bh_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('bh_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_eps():
        """ gravitational softening length, spline soft. (code length, 1.) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('eps', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_eps():
        """ gravitational softening length, spline soft. (code length, 1.) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('eps', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    def set_eps2(self,eps2):
        return self.set_eps(eps2**0.5)
    def get_eps2(self):
        eps,err=self.get_eps()
        return eps**2, err

    @legacy_function
    def set_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .01) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('gdgtol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .01) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('gdgtol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_nn_tol():
        """ fractional tolerance in nn_target  (0.1) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('nn_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_nn_tol():
        """ fractional tolerance in nn_target  (0.1) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('nn_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_epsgas():
        """ gas grav smoothing eps"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('epsgas', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_epsgas():
        """ gas grav smoothing eps"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('epsgas', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_gamma():
        """ gas polytropic index (1.666667) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gamma', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_gamma():
        """ gas polytropic index (1.666667) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gamma', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('alpha', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('alpha', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_beta():
        """ SPH artificial viscosity beta parameter (2*alpha=1.0) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('beta', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_beta():
        """ SPH artificial viscosity beta parameter (2*alpha=1.0) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('beta', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_epssph():
        """ SPH artificial viscosity safety against divergence (0.01) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('epssph', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_epssph():
        """ SPH artificial viscosity safety against divergence (0.01) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('epssph', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('courant', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('courant', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_removgas():
        """ minimum gas particle mass (fraction of initial (average) mass) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('removgas', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_removgas():
        """ minimum gas particle mass (fraction of initial (average) mass) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('removgas', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_consthsm():
        """ SPH smoothing length if constant"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('consthsm', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_consthsm():
        """ SPH smoothing length if constant"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('consthsm', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_nsmtol():
        """ tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmtol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_nsmtol():
        """ tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmtol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_graineff():
        """ FUV grain heating efficiency parameter (unitless, 0.05) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('graineff', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_graineff():
        """ FUV grain heating efficiency parameter (unitless, 0.05) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('graineff', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_crionrate():
        """ primary cosmic ray ionization rate (in units of 1.8e-17 sec^-1, 1.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('crionrate', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_crionrate():
        """ primary cosmic ray ionization rate (in units of 1.e-17 sec^-1, 3.6) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('crionrate', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_heat_par1():
        """ additional heating 1 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par1', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_heat_par1():
        """ additional heating 1 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par1', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_heat_par2():
        """ additional heating 2 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_heat_par2():
        """ additional heating 2 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_cool_par():
        """ additional cooling (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('cool_par', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_cool_par():
        """ additional cooling (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('cool_par', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_optdepth():
        """ 1/(mean free path) for UV photons (code length **-1, 0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('optdepth', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_optdepth():
        """ 1/(mean free path) for UV photons (code length **-1, 0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('optdepth', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_tcollfac():
        """ star formation delay parameter (unitless, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tcollfac', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_tcollfac():
        """ star formation delay parameter (unitless, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tcollfac', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_masscrit():
        """ star formation cloud reference mass (Msun, 1.e5) """
        function = LegacyFunctionSpecification()  
        function.addParameter('masscrit', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_masscrit():
        """ star formation cloud reference mass (Msun, 1.e5) """
        function = LegacyFunctionSpecification()  
        function.addParameter('masscrit', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_sfeff():
        """ gas particle mass fraction converted to stars (0.125) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfeff', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sfeff():
        """ gas particle mass fraction converted to stars (0.125) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfeff', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_tbubble():
        """ Supernova activity time, (Myr, 3.e7) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tbubble', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_tbubble():
        """ Supernova activity time, (Myr, 3.e7) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tbubble', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_sne_eff():
        """ Supernova feedback coupling efficiency, (0.1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sne_eff', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sne_eff():
        """ Supernova feedback coupling efficiency, (0.1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sne_eff', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_tsnbeg():
        """ Supernova feedback start time, (Myr, 3.e6) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tsnbeg', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_tsnbeg():
        """ Supernova feedback start time, (Myr, 3.e6) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tsnbeg', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_rhomax():
        """ Maximum permissible density (code dens units, 100) """
        function = LegacyFunctionSpecification()  
        function.addParameter('rhomax', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_rhomax():
        """ Maximum permissible density (code dens units, 100) """
        function = LegacyFunctionSpecification()  
        function.addParameter('rhomax', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function;


# character
    @legacy_function
    def set_halofile():
        """ halo model file (none) """
        function = LegacyFunctionSpecification()  
        function.addParameter('halofile', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_halofile():
        """ halo model file (none) """
        function = LegacyFunctionSpecification()  
        function.addParameter('halofile', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_feedback():
        """ feedback model (fuv, pres, kine, solo, solh) """
        function = LegacyFunctionSpecification()  
        function.addParameter('feedback', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_feedback():
        """ feedback model (fuv, pres, kine, solo, solh) """
        function = LegacyFunctionSpecification()  
        function.addParameter('feedback', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_sfmode():
        """ star formation model (gerritsen, nieuw) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfmode', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sfmode():
        """ star formation model (gerritsen, nieuw) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfmode', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_hupdatemethod():
        """ SPH smoothing length criterion (at the moment always 'mass')  """
        function = LegacyFunctionSpecification()  
        function.addParameter('hupdatemethod', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_hupdatemethod():
        """ SPH smoothing length criterion (at the moment always 'mass')  """
        function = LegacyFunctionSpecification()  
        function.addParameter('hupdatemethod', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_sph_visc():
        """ SPH viscosity (sph,sphv, bulk)  
        not all may work """
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_visc', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sph_visc():
        """ SPH viscosity (sph,sphv, bulk)  
        not all may work """
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_visc', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_time_step():
        """
        Set the model timestep.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('time_step', dtype='float64', direction=function.IN,
            description = "The current model timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time step was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function
    
    @legacy_function
    def set_fi_data_directory():
        """
        Update the path to the Fi database.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('fi_data_directory', dtype='string', direction=function.IN,
            description = "Name of the Fi data directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function
    
    @legacy_function
    def get_fi_data_directory():
        """
        Retrieve the path to the Fi database currently used.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('fi_data_directory', dtype='string', direction=function.OUT,
            description = "Name of the Fi data directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Value was retrieved
        -1 - ERROR
            Could not retrieve value
        """
        return function
    
    def get_periodic_boundaries_flag(self):
        return self.mode == self.MODE_PERIODIC_BOUNDARIES, 0
        
    @legacy_function
    def get_number_of_sph_particles_removed():
        """
        Return the number of particles deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32',
                              direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_id_of_removed_sph_particle():
        """
        Return the id of the nth particle deleted during the last evolve.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_delete', dtype='int32',
                              direction=function.IN, 
                 description = 'index in the deleted particles list (zero based)')
        function.addParameter('index_of_particle', dtype='int32',
                              direction=function.OUT)
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
        
    
class GlFiInterface(FiInterface):
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'fi_worker_gl'
        elif mode == self.MODE_PERIODIC_BOUNDARIES:
            return 'fi_worker_periodic_gl'
        else:
            return 'fi_worker_gl'

    def __init__(self,mode=FiInterface.MODE_NORMAL, **options):
        CodeInterface.__init__(self,name_of_the_worker = self.name_of_the_worker(mode), **options)
        
    @legacy_function
    def viewer():
        function = LegacyFunctionSpecification()  
        return function
    
    def start_viewer(self):
        self.viewer()
        
class FiDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__

class Fi(GravitationalDynamics):
    
    __doc__ = FiDoc()
    
    def __init__(self, convert_nbody = None, mode = 'normal', use_gl = False, **options):
        if(use_gl):
            legacy_interface = GlFiInterface(mode = mode, **options)
        else:
            legacy_interface = FiInterface(mode = mode, **options)            
        self.legacy_doc = legacy_interface.__doc__
        
        #if convert_nbody is None:
        #    convert_nbody=nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)

        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )     
    
    def initialize_code(self):
        result = self.overridden().initialize_code()
        
        self.parameters.fi_data_directory = self.legacy_interface.get_data_directory()+'/'
        
        if not self.unit_converter is None:
            value=self.unit_converter.to_si(nbody_system.length)
            self.parameters._original.code_length_unit = value
            
            value=self.unit_converter.to_si(nbody_system.mass)
            self.parameters._original.code_mass_unit = value
        
        return result

    def define_properties(self, object):
        GravitationalDynamics.define_properties(self, object)
        object.add_property("get_thermal_energy")
        object.add_property("get_total_energy")
    
    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        object.add_method('END', 'initialize_code')

        object.add_method('EDIT', 'new_dm_particle')
        object.add_method('UPDATE', 'new_dm_particle')
        object.add_transition('RUN', 'UPDATE', 'new_dm_particle', False)
        object.add_method('EDIT', 'new_sph_particle')
        object.add_method('UPDATE', 'new_sph_particle')
        object.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        object.add_method('EDIT', 'new_star_particle')
        object.add_method('UPDATE', 'new_star_particle')
        object.add_transition('RUN', 'UPDATE', 'new_star_particle', False)
        object.add_method('RUN', 'get_velocity')
        object.add_method('RUN', 'get_acceleration')
        object.add_method('RUN', 'get_internal_energy')
        object.add_method('RUN', 'get_dinternal_energy_dt')
        object.add_method('RUN', 'get_smoothing_length')
        object.add_method('RUN', 'get_density')
        object.add_method('RUN', 'get_star_tform')
        object.add_method('RUN', 'get_state_sph')
        object.add_method('RUN', 'get_state_star')
        
        object.add_method('RUN', 'get_kinetic_energy')
        object.add_method('RUN', 'get_potential_energy')
        object.add_method('RUN', 'get_thermal_energy')
        object.add_method('RUN', 'get_total_energy')
        object.add_method('RUN', 'get_total_radius')
        object.add_method('RUN', 'get_center_of_mass_position')
        object.add_method('RUN', 'get_center_of_mass_velocity')
        object.add_method('RUN', 'get_total_mass')
        object.add_method('RUN', 'get_time')
        object.add_method('EDIT', 'get_time')
        object.add_method('INITIALIZED', 'get_time')
        object.add_method('PARAMETER_CHANGE_A', 'get_time')
        object.add_method('PARAMETER_CHANGE_B', 'get_time')
        object.add_method('RUN', 'get_hydro_state_at_point')
    
# this should be checked!
        object.add_method('EDIT', 'get_gravity_at_point')
        object.add_method('EDIT', 'get_potential_at_point')
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2", 
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        
        object.add_method_parameter(
            "get_dtime", 
            "set_dtime",
            "timestep", 
            "timestep for system", 
            default_value = 1.0 | nbody_system.time
        ) 
        
        
        object.add_boolean_parameter(
            "get_radiate",
            "set_radiate",
            "radiation_flag",
            "Radiation flag. True means: radiation (i.e. radiative cooling/heating) is included. "
                "False means: no radiation, and implies no star formation.",
            False
        )
        
        object.add_boolean_parameter(
            "get_starform",
            "set_starform",
            "star_formation_flag",
            "Star-formation flag. True means: star formation is included. "
                "False means: no star formation included.",
            False
        )
        
        object.add_boolean_parameter(
            "get_use_hydro",
            "set_use_hydro",
            "use_hydro_flag",
            "Hydrodynamics flag. True means: SPH hydro included, False means: gravity only.",
            True
        )
        
        object.add_boolean_parameter(
            "get_sqrttstp",
            "set_sqrttstp",
            "square_root_timestep_flag",
            "Square-root-timestep flag. True means: use sqrt(eps/acc) timestep criterion.",
            False
        )
        
        object.add_boolean_parameter(
            "get_acc_tstp",
            "set_acc_tstp",
            "acc_timestep_flag",
            "Acceleration-timestep flag. True means: use vref/acc timestep criterion.",
            True
        )
        
        object.add_boolean_parameter(
            "get_freetstp",
            "set_freetstp",
            "freeform_timestep_flag",
            "Freeform-timestep flag. True means: use freeform (v/freev)**freevexp * (a/freea)**freeaexp timestep criterion.",
            False
        )
        
        object.add_boolean_parameter(
            "get_usequad",
            "set_usequad",
            "quadrupole_moments_flag",
            "Quadrupole-moments flag. True means: calculate and use quadrupole cell moments.",
            False
        )
        
        object.add_boolean_parameter(
            "get_directsum",
            "set_directsum",
            "direct_sum_flag",
            "Direct-summation flag. True means: direct N**2 gravity summation.",
            False
        )
        
        object.add_boolean_parameter(
            "get_selfgrav",
            "set_selfgrav",
            "self_gravity_flag",
            "Self-gravity flag. False means: self-gravity is not used, only external potentials.",
            True
        )
        
        object.add_boolean_parameter(
            "get_fixthalo",
            "set_fixthalo",
            "fixed_halo_flag",
            "Fixed-halo flag. True means: use fixed (spherical) potential.",
            False
        )
        
        object.add_boolean_parameter(
            "get_adaptive_eps",
            "set_adaptive_eps",
            "adaptive_smoothing_flag",
            "Adaptive-smoothing flag. True means: use of adaptive gravity smoothing for all particles.",
            False
        )
        
        object.add_boolean_parameter(
            "get_gdgop",
            "set_gdgop",
            "gadget_cell_opening_flag",
            "Gadget-cell-opening flag. True means: use of Gadget cell opening criterion.",
            True
        )
        
        object.add_boolean_parameter(
            "get_smoothinput",
            "set_smoothinput",
            "smooth_input_flag",
            "Smooth-input flag. True means: smooth input SPH properties.",
            False
        )
        
        object.add_boolean_parameter(
            "get_consph",
            "set_consph",
            "conservative_sph_flag",
            "Conservative-SPH flag. True means: use Springel & Hernquist conservative SPH form (currently the only option).",
            True
        )
        
        object.add_boolean_parameter(
            "get_sphinit",
            "set_sphinit",
            "sph_dens_init_flag",
            "SPH-density-init flag. True means: initialize sph density and h_smooth (most probably useless for AMUSE interface).",
            True
        )
        
        object.add_boolean_parameter(
            "get_uentropy",
            "set_uentropy",
            "integrate_entropy_flag",
            "Integrate-entropy flag. True means: integrate entropy, else: internal energy.",
            True
        )
        
        object.add_boolean_parameter(
            "get_isotherm",
            "set_isotherm",
            "isothermal_flag",
            "Isothermal flag. True means: isothermal gas (requires integrate_entropy_flag == False).",
            False
        )
        
        object.add_boolean_parameter(
            "get_eps_is_h",
            "set_eps_is_h",
            "eps_is_h_flag",
            "Eps-is-h flag. True means: set gas particles gravitational epsilon to h (SPH smoothing length).",
            True
        )
        
        
        object.add_method_parameter(
            "get_firstsnap", 
            "set_firstsnap",
            "first_snapshot", 
            "The number of the first snapshot.", 
            default_value = 0
        )
        
        object.add_method_parameter(
            "get_stepout", 
            "set_stepout",
            "output_interval", 
            "The number of steps between output.", 
            default_value = 5
        )
        
        object.add_method_parameter(
            "get_steplog", 
            "set_steplog",
            "log_interval", 
            "The number of steps between logs.", 
            default_value = 5
        )
        
        object.add_method_parameter(
            "get_max_tbin", 
            "set_max_tbin",
            "maximum_time_bin", 
            "The maximum time bin (dtime*2**-max_tbin=minimum time step).", 
            default_value = 4096
        )
        
        object.add_method_parameter(
            "get_minppbin", 
            "set_minppbin",
            "minimum_part_per_bin", 
            "The minimum number of particles per time bin.", 
            default_value = 1
        )
        
        object.add_method_parameter(
            "get_targetnn", 
            "set_targetnn",
            "targetnn", 
            "The target number of neighbour particles for variable gravitational eps.", 
            default_value = 32
        )
        
        object.add_method_parameter(
            "get_verbosity", 
            "set_verbosity",
            "verbosity", 
            "The level of terminal output (0=minimum).", 
            default_value = 0
        )
        
        object.add_method_parameter(
            "get_nsmooth", 
            "set_nsmooth",
            "n_smooth", 
            "The target number of SPH neighbours.", 
            default_value = 64
        )
        
        
        object.add_method_parameter(
            "get_pboxsize", 
            "set_pboxsize",
            "periodic_box_size", 
            "The size of simulation domain box (particles outside get deleted).", 
            default_value = 300.0 | nbody_system.length
        )
        
        object.add_method_parameter(
            "get_unitm_in_msun", 
            "set_unitm_in_msun",
            "code_mass_unit", 
            "The code mass unit (in Msun, 10^9 standard).", 
            default_value = 1.0e9 | units.MSun
        )
        
        object.add_method_parameter(
            "get_unitl_in_kpc", 
            "set_unitl_in_kpc",
            "code_length_unit", 
            "The code length unit (in kpc, 1 standard).", 
            default_value = 1.0 | units.kpc
        )
        
        object.add_method_parameter(
            "get_tstepcrit", 
            "set_tstepcrit",
            "sqrt_timestep_crit_constant", 
            "Square-root-timestep criterion constant (unitless,standard=1.).", 
            default_value = 1.0
        )
        
        object.add_method_parameter(
            "get_tstpcr2", 
            "set_tstpcr2",
            "acc_timestep_crit_constant", 
            "Acceleration-timestep criterion constant (unitless,standard=0.25).", 
            default_value = 0.25
        )
        
        object.add_method_parameter(
            "get_freev", 
            "set_freev",
            "free_timestep_crit_constant_v", 
            "Freeform-timestep criterion constant v.", 
            default_value = 0.5
        )
        
        object.add_method_parameter(
            "get_freea", 
            "set_freea",
            "free_timestep_crit_constant_a", 
            "Freeform-timestep criterion constant a.", 
            default_value = 0.35
        )
        
        object.add_method_parameter(
            "get_freevexp", 
            "set_freevexp",
            "free_timestep_crit_constant_vexp", 
            "Freeform-timestep criterion constant v_exp.", 
            default_value = 0.0
        )
        
        object.add_method_parameter(
            "get_freeaexp", 
            "set_freeaexp",
            "free_timestep_crit_constant_aexp", 
            "Freeform-timestep criterion constant a_exp.", 
            default_value = -1.0
        )
        
        object.add_method_parameter(
            "get_bh_tol", 
            "set_bh_tol",
            "opening_angle", 
            "Opening angle, theta, for building the tree: between 0 and 1 (unitless, 0.5).", 
            default_value = 0.5
        )
        
        object.add_method_parameter(
            "get_gdgtol", 
            "set_gdgtol",
            "gadget_cell_opening_constant", 
            "Gadget-cell-openings criterion parameter  (unitless, .01)", 
            default_value = 0.01
        )
        
        object.add_method_parameter(
            "get_nn_tol", 
            "set_nn_tol",
            "nn_tol", 
            "The fractional tolerance in nn_target (0.1).", 
            default_value = 0.1
        )
        
        object.add_method_parameter(
            "get_epsgas", 
            "set_epsgas",
            "gas_epsilon", 
            "The gas gravitational smoothing epsilon.", 
            default_value = 0.005 | nbody_system.length
        )
        
        object.add_method_parameter(
            "get_gamma", 
            "set_gamma",
            "gamma", 
            "gas polytropic index (1.6666667)", 
            default_value = 1.6666667
        )
        
        object.add_method_parameter(
            "get_alpha", 
            "set_alpha",
            "artificial_viscosity_alpha", 
            "SPH artificial viscosity alpha parameter (0.5)", 
            default_value = 0.5
        )
        
        object.add_method_parameter(
            "get_beta", 
            "set_beta",
            "beta", 
            "SPH artificial viscosity beta parameter (2*alpha=1.0)", 
            default_value = 1.0
        )
        
        object.add_method_parameter(
            "get_epssph", 
            "set_epssph",
            "sph_artificial_viscosity_eps", 
            "SPH artificial viscosity safety against divergence (0.01)", 
            default_value = 0.01
        )
        
        object.add_method_parameter(
            "get_courant", 
            "set_courant",
            "courant", 
            "SPH courant condition parameter (0.3)", 
            default_value = 0.3
        )
        
        object.add_method_parameter(
            "get_removgas", 
            "set_removgas",
            "min_gas_part_mass", 
            "minimum gas particle mass (fraction of initial (average) mass)", 
            default_value = 0.25
        )
        
        object.add_method_parameter(
            "get_consthsm", 
            "set_consthsm",
            "sph_h_const", 
            "SPH smoothing length if constant", 
            default_value = 0.2 | nbody_system.length
        )
        
        object.add_method_parameter(
            "get_nsmtol", 
            "set_nsmtol",
            "n_smooth_tol", 
            "fractional tolerance in number of SPH neighbours", 
            default_value = 0.1
        )
        
        object.add_method_parameter(
            "get_graineff", 
            "set_graineff",
            "grain_heat_eff", 
            "FUV grain heating efficiency parameter (unitless, 0.05)", 
            default_value = 0.05
        )
        
        object.add_method_parameter(
            "get_crionrate", 
            "set_crionrate",
            "zeta_cr_ion_rate", 
            "primary cosmic ray ionization rate, zeta (in units of 1.8e-17 sec^-1, 1.)", 
            default_value = 3.6 | 1.8e-17 * units.s**-1
        )
        
        object.add_method_parameter(
            "get_heat_par1", 
            "set_heat_par1",
            "heat_par1", 
            "additional heating 1 parameter (0.0)", 
            default_value = 0.0
        )
        
        object.add_method_parameter(
            "get_heat_par2", 
            "set_heat_par2",
            "heat_par2", 
            "additional heating 2 parameter (0.0)", 
            default_value = 0.0
        )
        
        object.add_method_parameter(
            "get_cool_par", 
            "set_cool_par",
            "cool_par", 
            "additional cooling parameter (1.0)", 
            default_value = 1.0
        )
        
        object.add_method_parameter(
            "get_optdepth", 
            "set_optdepth",
            "optical_depth", 
            "1/(mean free path) for UV photons (code length **-1, 0.0)", 
            default_value = 0.0
        )
        
        object.add_method_parameter(
            "get_tcollfac", 
            "set_tcollfac",
            "star_form_delay_fac", 
            "star formation delay parameter (unitless, 1)", 
            default_value = 1.0
        )
        
        object.add_method_parameter(
            "get_masscrit", 
            "set_masscrit",
            "star_form_mass_crit", 
            "star formation cloud reference mass (Msun, 1.e5)", 
            default_value = 1.0e5 | units.MSun
        )
        
        object.add_method_parameter(
            "get_sfeff", 
            "set_sfeff",
            "star_form_eff", 
            "gas particle mass fraction converted to stars (0.125)", 
            default_value = 0.25
        )
        
        object.add_method_parameter(
            "get_tbubble", 
            "set_tbubble",
            "supernova_duration", 
            "Supernova activity time, (Myr, 3.e7)", 
            default_value = 3.0e7 | units.Myr
        )
        
        object.add_method_parameter(
            "get_sne_eff", 
            "set_sne_eff",
            "supernova_eff", 
            "Supernova feedback coupling efficiency, (0.0)", 
            default_value = 0.0
        )
        
        object.add_method_parameter(
            "get_tsnbeg", 
            "set_tsnbeg",
            "t_supernova_start", 
            "Supernova feedback start time, (Myr, 3.e6)", 
            default_value = 3.0e6 | units.Myr
        )
        
        object.add_method_parameter(
            "get_rhomax", 
            "set_rhomax",
            "max_density", 
            "Maximum permissible density (code density units, 100)", 
            default_value = 100.0
        )
        
        object.add_method_parameter(
            "get_halofile", 
            "set_halofile",
            "halofile", 
            "Path to initial halo model file, relative to the Fi data directory (none)", 
            default_value = "none"
        )
        
        object.add_method_parameter(
            "get_feedback", 
            "set_feedback",
            "feedback", 
            "feedback model (fuv, pres, kine, solo, solh)", 
            default_value = "fuv"
        )
        
        object.add_method_parameter(
            "get_sfmode", 
            "set_sfmode",
            "star_formation_mode", 
            "star formation model (gerritsen, nieuw)", 
            default_value = "gerritsen"
        )
        
        object.add_method_parameter(
            "get_hupdatemethod", 
            "set_hupdatemethod",
            "h_update_method", 
            "SPH smoothing length criterion (at the moment always 'mass')", 
            default_value = "mass"
        )
        
        object.add_method_parameter(
            "get_sph_visc", 
            "set_sph_visc",
            "sph_viscosity", 
            "SPH viscosity (sph,sphv, bulk). Note: not all may work.", 
            default_value = "sph"
        )
        
        object.add_method_parameter(
            "get_fi_data_directory", 
            "set_fi_data_directory",
            "fi_data_directory", 
            "Name of the Fi data directory", 
            default_value = ""
        )

        object.add_boolean_parameter(
            "get_periodic_boundaries_flag",
            None,
            "periodic_boundaries_flag",
            "Periodic boundaries flag. True means: use periodic boundary conditions (read-only)",
            False
        )
        
        self.stopping_conditions.define_parameters(object)        
    
    def define_particle_sets(self, object):
        object.define_super_set('particles', ['dm_particles','gas_particles','star_particles'], 
            index_to_default_set = 0)
        
        object.define_set('dm_particles', 'id')
        object.set_new('dm_particles', 'new_dm_particle')
        object.set_delete('dm_particles', 'delete_particle')
        object.add_setter('dm_particles', 'set_state')
        object.add_getter('dm_particles', 'get_state')
        object.add_setter('dm_particles', 'set_mass')
        object.add_getter('dm_particles', 'get_mass', names = ('mass',))
        object.add_setter('dm_particles', 'set_position')
        object.add_getter('dm_particles', 'get_position')
        object.add_setter('dm_particles', 'set_radius')
        object.add_getter('dm_particles', 'get_radius')
        object.add_setter('dm_particles', 'set_velocity')
        object.add_getter('dm_particles', 'get_velocity')
        
        object.define_set('gas_particles', 'id')
        object.set_new('gas_particles', 'new_sph_particle')
        object.set_delete('gas_particles', 'delete_particle')
        object.add_setter('gas_particles', 'set_state_sph')
        object.add_getter('gas_particles', 'get_state_sph')
        object.add_setter('gas_particles', 'set_mass')
        object.add_getter('gas_particles', 'get_mass', names = ('mass',))
        object.add_getter('gas_particles', 'get_radius')
        object.add_setter('gas_particles', 'set_radius')
        object.add_setter('gas_particles', 'set_position')
        object.add_getter('gas_particles', 'get_position')
        object.add_setter('gas_particles', 'set_velocity')
        object.add_getter('gas_particles', 'get_velocity')
        object.add_setter('gas_particles', 'set_internal_energy')
        object.add_getter('gas_particles', 'get_internal_energy')
        object.add_getter('gas_particles', 'get_dinternal_energy_dt')
        object.add_setter('gas_particles', 'set_smoothing_length')
        object.add_getter('gas_particles', 'get_smoothing_length')
        object.add_getter('gas_particles', 'get_density')
        
        object.define_set('star_particles', 'id')
        object.set_new('star_particles', 'new_star_particle')
        object.set_delete('star_particles', 'delete_particle')
        object.add_setter('star_particles', 'set_state_star')
        object.add_getter('star_particles', 'get_state_star')
        object.add_setter('star_particles', 'set_mass')
        object.add_getter('star_particles', 'get_mass', names = ('mass',))
        object.add_setter('star_particles', 'set_position')
        object.add_getter('star_particles', 'get_position')
        object.add_setter('star_particles', 'set_radius')
        object.add_getter('star_particles', 'get_radius')
        object.add_setter('star_particles', 'set_velocity')
        object.add_getter('star_particles', 'get_velocity')
        object.add_setter('star_particles', 'set_star_tform')
        object.add_getter('star_particles', 'get_star_tform')
        
        self.stopping_conditions.define_particle_set(object, 'particles')

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        object.add_method(
            "set_velocity",
            (
                object.INDEX,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_velocity",
            (
                object.INDEX,
            ),
            (
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "new_dm_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "new_sph_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_state_sph",
            (
                object.INDEX,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state_sph",
            (
                object.INDEX,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "set_internal_energy",
            (
                object.INDEX,
                nbody_system.specific_energy,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_internal_energy",
            (
                object.INDEX,
            ),
            (
                nbody_system.specific_energy,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_dinternal_energy_dt",
            (
                object.INDEX,
            ),
            (
                nbody_system.specific_energy/nbody_system.time,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_smoothing_length",
            (object.INDEX, nbody_system.length),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_smoothing_length",
            (object.INDEX,),
            (nbody_system.length, object.ERROR_CODE)
        )
        object.add_method(
            "get_density",
            (object.INDEX,),
            (nbody_system.density, object.ERROR_CODE)
        )
        
        object.add_method(
            "new_star_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.time,
                nbody_system.length,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_state_star",
            (
                object.INDEX,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.time,
                nbody_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state_star",
            (
                object.INDEX,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.time,
                nbody_system.length,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "set_star_tform",
            (
                object.INDEX,
                nbody_system.time,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_star_tform",
            (
                object.INDEX,
            ),
            (
                nbody_system.time,
                object.ERROR_CODE
            )
        )
        
        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration, nbody_system.acceleration, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_hydro_state_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length,
                nbody_system.speed, nbody_system.speed, nbody_system.speed),
            (nbody_system.density, nbody_system.momentum_density, nbody_system.momentum_density, 
                nbody_system.momentum_density, nbody_system.energy_density, object.ERROR_CODE)
        )
        
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_dtime",
            (),
            (nbody_system.time, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_dtime",
            (nbody_system.time, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_firstsnap",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_firstsnap",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_stepout",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_stepout",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_steplog",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_steplog",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_max_tbin",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_max_tbin",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_minppbin",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_minppbin",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_targetnn",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_targetnn",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_verbosity",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_verbosity",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_nsmooth",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_nsmooth",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_pboxsize",
            (),
            (nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_pboxsize",
            (nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_unitm_in_msun",
            (),
            (units.MSun, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_unitm_in_msun",
            (units.MSun, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_unitl_in_kpc",
            (),
            (units.kpc, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_unitl_in_kpc",
            (units.kpc, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tstepcrit",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_tstepcrit",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tstpcr2",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_tstpcr2",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_freev",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_freev",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_freea",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_freea",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_freevexp",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_freevexp",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_freeaexp",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_freeaexp",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_bh_tol",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_bh_tol",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_gdgtol",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_gdgtol",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_nn_tol",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_nn_tol",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_epsgas",
            (),
            (nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_epsgas",
            (nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_gamma",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_gamma",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_alpha",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_alpha",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_beta",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_beta",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_epssph",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_epssph",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_courant",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_courant",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_removgas",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_removgas",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_consthsm",
            (),
            (nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_consthsm",
            (nbody_system.length, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_nsmtol",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_nsmtol",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_graineff",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_graineff",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_crionrate",
            (),
            (1.8e-17 * units.s**-1, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_crionrate",
            (1.8e-17 * units.s**-1, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_heat_par1",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_heat_par1",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_heat_par2",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_heat_par2",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_cool_par",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_cool_par",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_optdepth",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_optdepth",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tcollfac",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_tcollfac",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_masscrit",
            (),
            (units.MSun, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_masscrit",
            (units.MSun, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_sfeff",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_sfeff",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tbubble",
            (),
            (units.Myr, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_tbubble",
            (units.Myr, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_sne_eff",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_sne_eff",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_tsnbeg",
            (),
            (units.Myr, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_tsnbeg",
            (units.Myr, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_rhomax",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_rhomax",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_halofile",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_halofile",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_feedback",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_feedback",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_sfmode",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_sfmode",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_hupdatemethod",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_hupdatemethod",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_sph_visc",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_sph_visc",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_fi_data_directory",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_fi_data_directory",
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_thermal_energy",
            (),
            (nbody_system.energy, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_total_energy",
            (),
            (nbody_system.energy, object.ERROR_CODE,)
        )
        
        object.add_method(
            "get_number_of_sph_particles_removed",
            (
            ),
            (
                object.NO_UNIT,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_id_of_removed_sph_particle",
            (
                object.NO_UNIT,
            ),
            (
                object.INDEX,
                object.ERROR_CODE
            )
        )

        self.stopping_conditions.define_methods(object)       
        
    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.datamodel.incode_storage module, as
        it uses a lot of internal methods and info!
        """
        number_of_removed_particles = self.get_number_of_sph_particles_removed()
        
        if number_of_removed_particles == 0:
            return
        
        indices_in_update_list = range(number_of_removed_particles)
        indices_to_remove = self.get_id_of_removed_sph_particle(indices_in_update_list)
        
        incode_storage = self.gas_particles._private.attribute_storage
        incode_storage._remove_indices(indices_to_remove)
