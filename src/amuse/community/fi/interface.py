import os
import numpy
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode
from amuse.community import *
from amuse.support.options import option

from amuse.units import generic_unit_system 
from amuse.community.interface.common import CommonCode


class FiInterface(
        CodeInterface,
        GravitationalDynamicsInterface,
        LiteratureReferencesMixIn,
        StoppingConditionInterface,
        GravityFieldInterface,
        CodeWithDataDirectories):   
    """
    FI is a parallel TreeSPH code for galaxy simulations. Extensively 
    rewritten, extended and parallelized it is a development from code from 
    Jeroen Gerritsen and Roelof Bottema, which itself goes back to Treesph. 
    
    Note that some features are not working atm. These may be fixed in the 
    future. (and I will think of a better name)
    
    The relevant references are:
        .. [#] ** Hernquist & Katz, 1989, ApJS 70, 419 [1989ApJS...70..419H]
        .. [#] ** Pelupessy, van der Werf & Icke, 2004, A&A 422, 55 [2004A&A...422...55P]
        .. [#] Pelupessy, PhD thesis 2005, Leiden Observatory [2005PhDT........17P]
        .. [#] Gerritsen & Icke, 1997, A&A 325, 972 [1997A&A...325..972G]
    """
    get_total_radius=None
    get_total_mass=None
    get_center_of_mass_position=None
    get_center_of_mass_velocity=None
    set_acceleration=None
    get_acceleration=None
    
    use_modules=['StoppingConditions','AmuseInterface']
    
    MODE_NORMAL = 'normal'
    MODE_NORMAL_OPENMP = 'openmp'
    MODE_PERIODIC_BOUNDARIES = 'periodic'
    
    def __init__(self, mode = MODE_NORMAL,  **options):
        self.mode = mode
        
        if "number_of_workers" in options and options["number_of_workers"]!=1:
          raise Exception("Fi expects number_of_workers to be 1,\nfor multiple processors use mode='openmp' and set the OMP_NUM_THREADS environment variable") 
        
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(mode), **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
                     
    
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'fi_worker'
        elif mode == self.MODE_NORMAL_OPENMP:
            return 'fi_worker_mp'
        elif mode == self.MODE_PERIODIC_BOUNDARIES:
            return 'fi_worker_periodic'
        else:
            return 'fi_worker'
    
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
    def get_pressure():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='int32', direction=function.IN)
        function.addParameter('pressure', dtype='float64', direction=function.OUT)
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
    def synchronize_model():
        """ synchronize the model """
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function
    def trigger_partremoval():
        function = LegacyFunctionSpecification()  
        return function


# setting/ getting parameters
# logicals
    @legacy_function   
    def set_use_hydro():
        """ set_use_hydro([0,1]): SPH hydro if 1, gravity only if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('use_hydro_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_use_hydro():
        """ get_use_hydro(): SPH hydro if 1, gravity only if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('use_hydro_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_radiate():
        """ set_radiate([0,1]): rad cooling/heating if 1, not if 0
             radiate false implies starform false """    
        function = LegacyFunctionSpecification()  
        function.addParameter('radiation_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_radiate():
        """ get_radiate(): rad cooling/heating if 1, not if 0
             radiate false implies starform false """        
        function = LegacyFunctionSpecification()  
        function.addParameter('radiation_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_starform():
        """ set_starform([0,1]): star formation if 1, not if 0 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('star_formation_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_starform():
        """ get_starform(): star formation if 1, not if 0 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('star_formation_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_cosmo():
        """ set_cosmo([0,1]): not functional at the moment """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_cosmo():
        """ get_cosmo(): not functional at the moment """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_sqrttstp():
        """ set_sqrttstp([0,1]): use sqrt(eps/acc) timestep crit if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('square_root_timestep_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_sqrttstp():
        """ get_sqrttstp(): use sqrt(eps/acc) timestep crit if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('square_root_timestep_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_acc_tstp():
        """ set_acc_tstp([0,1]): use vref/acc timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('acc_timestep_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_acc_tstp():
        """ get_acc_tstp(): use vref/acc timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('acc_timestep_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_freetstp():
        """ set_freetstp([0,1]): use freeform (v/freev)**freevexp * (a/freea)**freeaexp timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('freeform_timestep_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_freetstp():
        """ get_freetstp(): use freeform timestep crit if 1"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('freeform_timestep_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_usequad():
        """ set_usequad([0,1]): calc. and use quadrupole cell moments if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('quadrupole_moments_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_usequad():
        """ get_usequad(): calc. and use quadrupole cell moments if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('quadrupole_moments_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_directsum():
        """ set_directsum([0,1]): direct N**2 grav sum if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('direct_sum_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_directsum():
        """ get_directsum(): direct N**2 grav sum if 1"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('direct_sum_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_selfgrav():
        """ set_selfgrav([0,1]): calculate self-gravity if 1
          if set to 0, self gravity is not used, only external potentials"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('self_gravity_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_selfgrav():
        """ get_selfgrav(): calculate self-gravity if 1 """
        function = LegacyFunctionSpecification()  
        function.addParameter('self_gravity_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_fixthalo():
        """ set_fixthalo([0,1]): use fixed (spherical) potential if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('fixed_halo_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_fixthalo():
        """ get_fixthalo(): use fixed (spherical) potential if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('fixed_halo_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_adaptive_eps():
        """ set_adaptive_eps([0,1]): use of adaptive grav smoothing for all part if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('adaptive_smoothing_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_adaptive_eps():
        """ get_adaptive_eps(): use of adaptive grav smoothing for all part if 1 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('adaptive_smoothing_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_gdgop():
        """ set_gdgop([0,1]): use of gadget cell opening criterion if 1 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_gdgop():
        """ get_gdgop(): use of gadget cell opening criterion if 1 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gadget_cell_opening_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_smoothinput():
        """ set_smoothinput([0,1]): smooth input SPH prop. if 1 
         (not working) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('smooth_input_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_smoothinput():
        """ get_smoothinput(): smooth input SPH prop. if 1 
         (not working) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('smooth_input_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_consph():
        """ set_consph([0,1]): use springel&Hernquist conservative SPH form. if 1 
          at the moment this is only option"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('conservative_sph_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_consph():
        """ get_consph(): use springel&Hernquist conservative SPH form. if 1 
          at the moment this is only option"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('conservative_sph_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_sphinit():
        """ set_sphinit([0,1]): initialize sph dens and hsmooth if 1 
         most probably useless for AMUSE interface"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_dens_init_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_sphinit():
        """ set_sphinit([0,1]): initialize sph dens and hsmooth if 1 
         most probably useless for AMUSE interface"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_dens_init_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_uentropy():
        """ set_uentropy([0,1]): integrate entropy if 1, internal energy if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('integrate_entropy_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_uentropy():
        """ get_uentropy(): integrate entropy if 1, internal energy if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('integrate_entropy_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_isotherm():
        """ set_isotherm([0,1]): isothermal gas if 1
          note that isotherm needs set_uentropy(0) (false)"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_isotherm():
        """ get_isotherm(): isothermal gas if 1
          note that isotherm needs set_uentropy(0) (false)"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('isothermal_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_eps_is_h():
        """ set_eps_is_h([0,1]): gas particles grav. eps to SPH h if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('eps_is_h_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_eps_is_h():
        """ get_eps_is_h(): gas particles grav. eps to SPH h if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('eps_is_h_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_balsara():
        """ set_balsara([0,1]): use Balsara viscosity limiter if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('balsara_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_balsara():
        """ set_balsara([0,1]): use Balsara viscosity limiter if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('balsara_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_mingaseps():
        """ set_mingaseps([0,1]): enforce minimum gas grav eps if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('mingaseps_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_mingaseps():
        """ set_mingaseps([0,1]): enforce minimum gas grav eps if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('mingaseps_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function



# integers
    @legacy_function
    def set_firstsnap():
        """ no. of first snapshot """
        function = LegacyFunctionSpecification()  
        function.addParameter('first_snapshot', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_firstsnap():
        """ no. of first snapshot """
        function = LegacyFunctionSpecification()  
        function.addParameter('first_snapshot', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_stepout():
        """ no. of steps between output """
        function = LegacyFunctionSpecification()  
        function.addParameter('output_interval', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_stepout():
        """ no. of steps between output """
        function = LegacyFunctionSpecification()  
        function.addParameter('output_interval', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_steplog():
        """ no. of steps between logs """
        function = LegacyFunctionSpecification()  
        function.addParameter('log_interval', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_steplog():
        """ no. of steps between logs """
        function = LegacyFunctionSpecification()  
        function.addParameter('log_interval', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_max_tbin():
        """ maximum time bin (dtime*2**-max_tbin=minimum time step)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('maximum_time_bin', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_max_tbin():
        """ maximum time bin (dtime*2**-max_tbin=minimum time step)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('maximum_time_bin', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_minppbin():
        """ target no. of particles per time bin"""
        function = LegacyFunctionSpecification()  
        function.addParameter('minimum_part_per_bin', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_minppbin():
        """ target no. of particles per time bin"""
        function = LegacyFunctionSpecification()  
        function.addParameter('minimum_part_per_bin', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_targetnn():
        """ target no. of neighbour particles for variable grav. eps"""
        function = LegacyFunctionSpecification()  
        function.addParameter('targetnn', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_targetnn():
        """ target no. of neighbour particles for variable grav. eps"""
        function = LegacyFunctionSpecification()  
        function.addParameter('targetnn', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_verbosity():
        """ level of terminal output (0=minimum)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('verbosity', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_verbosity():
        """ level of terminal output (0=minimum)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('verbosity', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmooth', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_nsmooth():
        """ target number of SPH neighbours"""
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmooth', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

# real
    @legacy_function
    def set_pboxsize():
        """ size of simulation domain box (particles outside get deleted)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('pboxsize', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_pboxsize():
        """ size of simulation domain box (particles outside get deleted)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('pboxsize', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_dtime():
        """ timestep (code units)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('dtime', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_dtime():
        """ timestep (code units)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('dtime', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_unitm_in_msun():
        """ code mass unit (in Msun, 10^9 standard) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitm_in_msun', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_unitm_in_msun():
        """ code mass unit (in Msun, 10^9 standard) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitm_in_msun', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_unitl_in_kpc():
        """ code length unit (in kpc, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitl_in_kpc', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_unitl_in_kpc():
        """ code length unit (in kpc, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('unitl_in_kpc', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function
    def set_tstepcrit():
        """ sqrttstp timestep constant (unitless,standard=1.) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstepcrit', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_tstepcrit():
        """ sqrttstp timestep constant (unitless,standard=1.) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstepcrit', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_tstpcr2():
        """ acc_tstp timestep constant (unitless,standard=0.25) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstpcr2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_tstpcr2():
        """ acc_tstp timestep constant (unitless,standard=0.25) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tstpcr2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_freev():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freev', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_freev():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freev', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_freea():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freea', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_freea():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freea', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_freevexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freevexp', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_freevexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freevexp', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_freeaexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freeaexp', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_freeaexp():
        """ freeform timestep constant """    
        function = LegacyFunctionSpecification()  
        function.addParameter('freeaexp', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('bh_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_bh_tol():
        """ Barnes Hut opening angle parameter (unitless, 0.5) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('bh_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_eps():
        """ gravitational softening length, spline soft. (code length, 1.) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('eps', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_eps():
        """ gravitational softening length, spline soft. (code length, 1.) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('eps', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

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
        return function
    @legacy_function   
    def get_gdgtol():
        """ Gadget cell openings criterion parameter  (unitless, .01) """    
        function = LegacyFunctionSpecification()  
        function.addParameter('gdgtol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_nn_tol():
        """ fractional tolerance in nn_target  (0.1) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('nn_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_nn_tol():
        """ fractional tolerance in nn_target  (0.1) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('nn_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_epsgas():
        """ gas grav smoothing eps"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('epsgas', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_epsgas():
        """ gas grav smoothing eps"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('epsgas', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_gamma():
        """ gas polytropic index (1.666667) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gamma', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_gamma():
        """ gas polytropic index (1.666667) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('gamma', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('alpha', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_alpha():
        """ SPH artificial viscosity alpha parameter (0.5) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('alpha', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_beta():
        """ SPH artificial viscosity beta parameter (2*alpha=1.0) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('beta', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_beta():
        """ SPH artificial viscosity beta parameter (2*alpha=1.0) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('beta', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_epssph():
        """ SPH artificial viscosity safety against divergence (0.01) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('epssph', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_epssph():
        """ SPH artificial viscosity safety against divergence (0.01) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('epssph', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('courant', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_courant():
        """ SPH courant condition parameter (0.3) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('courant', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_removgas():
        """ minimum gas particle mass (fraction of initial (average) mass) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('removgas', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_removgas():
        """ minimum gas particle mass (fraction of initial (average) mass) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('removgas', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_consthsm():
        """ SPH smoothing length if constant"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('consthsm', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_consthsm():
        """ SPH smoothing length if constant"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('consthsm', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_nsmtol():
        """ tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmtol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_nsmtol():
        """ tolerance in number of SPH neighbours """            
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmtol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_graineff():
        """ FUV grain heating efficiency parameter (unitless, 0.05) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('graineff', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_graineff():
        """ FUV grain heating efficiency parameter (unitless, 0.05) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('graineff', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_crionrate():
        """ primary cosmic ray ionization rate (in units of 1.8e-17 sec^-1, 1.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('crionrate', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_crionrate():
        """ primary cosmic ray ionization rate (in units of 1.e-17 sec^-1, 3.6) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('crionrate', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_heat_par1():
        """ additional heating 1 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par1', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_heat_par1():
        """ additional heating 1 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par1', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_heat_par2():
        """ additional heating 2 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_heat_par2():
        """ additional heating 2 (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('heat_par2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_cool_par():
        """ additional cooling (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('cool_par', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_cool_par():
        """ additional cooling (0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('cool_par', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_optdepth():
        """ 1/(mean free path) for UV photons (code length **-1, 0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('optdepth', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_optdepth():
        """ 1/(mean free path) for UV photons (code length **-1, 0.0)"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('optdepth', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_tcollfac():
        """ star formation delay parameter (unitless, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tcollfac', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_tcollfac():
        """ star formation delay parameter (unitless, 1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tcollfac', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_masscrit():
        """ star formation cloud reference mass (Msun, 1.e5) """
        function = LegacyFunctionSpecification()  
        function.addParameter('masscrit', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_masscrit():
        """ star formation cloud reference mass (Msun, 1.e5) """
        function = LegacyFunctionSpecification()  
        function.addParameter('masscrit', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_sfeff():
        """ gas particle mass fraction converted to stars (0.125) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfeff', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_sfeff():
        """ gas particle mass fraction converted to stars (0.125) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfeff', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_tbubble():
        """ Supernova activity time, (Myr, 3.e7) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tbubble', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_tbubble():
        """ Supernova activity time, (Myr, 3.e7) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tbubble', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_sne_eff():
        """ Supernova feedback coupling efficiency, (0.1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sne_eff', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_sne_eff():
        """ Supernova feedback coupling efficiency, (0.1) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sne_eff', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_tsnbeg():
        """ Supernova feedback start time, (Myr, 3.e6) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tsnbeg', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_tsnbeg():
        """ Supernova feedback start time, (Myr, 3.e6) """
        function = LegacyFunctionSpecification()  
        function.addParameter('tsnbeg', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_rhomax():
        """ Maximum density in case of star formation (force SF if exceeded, ignored if star formation is off) """
        function = LegacyFunctionSpecification()  
        function.addParameter('rhomax', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_rhomax():
        """ Maximum density in case of star formation (force SF if exceeded, ignored if star formation is off) """
        function = LegacyFunctionSpecification()  
        function.addParameter('rhomax', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function


# character
    @legacy_function
    def set_halofile():
        """ halo model file (none) """
        function = LegacyFunctionSpecification()  
        function.addParameter('halofile', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_halofile():
        """ halo model file (none) """
        function = LegacyFunctionSpecification()  
        function.addParameter('halofile', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_feedback():
        """ feedback model (fuv, pres, kine, solo, solh) """
        function = LegacyFunctionSpecification()  
        function.addParameter('feedback', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_feedback():
        """ feedback model (fuv, pres, kine, solo, solh) """
        function = LegacyFunctionSpecification()  
        function.addParameter('feedback', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_sfmode():
        """ star formation model (gerritsen, nieuw) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfmode', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_sfmode():
        """ star formation model (gerritsen, nieuw) """
        function = LegacyFunctionSpecification()  
        function.addParameter('sfmode', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_hupdatemethod():
        """ SPH smoothing length criterion (at the moment always 'mass')  """
        function = LegacyFunctionSpecification()  
        function.addParameter('hupdatemethod', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_hupdatemethod():
        """ SPH smoothing length criterion (at the moment always 'mass')  """
        function = LegacyFunctionSpecification()  
        function.addParameter('hupdatemethod', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_sph_visc():
        """ SPH viscosity (sph,sphv, bulk)  
        not all may work """
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_visc', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_sph_visc():
        """ SPH viscosity (sph,sphv, bulk)  
        not all may work """
        function = LegacyFunctionSpecification()  
        function.addParameter('sph_visc', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

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
        return self.mode == self.MODE_PERIODIC_BOUNDARIES
        
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
        self.mode=mode
        CodeInterface.__init__(self,name_of_the_worker = self.name_of_the_worker(mode), **options)

    @legacy_function   
    def get_image_target():
        """ target point of image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT, unit=nbody_system.length)
        function.addParameter('y', dtype='d', direction=function.OUT, unit=nbody_system.length)
        function.addParameter('z', dtype='d', direction=function.OUT, unit=nbody_system.length)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_image_target():
        """ target point of image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.IN, unit=nbody_system.length)
        function.addParameter('y', dtype='d', direction=function.IN, unit=nbody_system.length)
        function.addParameter('z', dtype='d', direction=function.IN, unit=nbody_system.length)
        function.result_type = 'i'
        return function

    @legacy_function   
    def get_viewpoint():
        """ camera position (for perspective proj) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT, unit=nbody_system.length)
        function.addParameter('y', dtype='d', direction=function.OUT, unit=nbody_system.length)
        function.addParameter('z', dtype='d', direction=function.OUT, unit=nbody_system.length)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_viewpoint():
        """ camera position (for perspective proj) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.IN, unit=nbody_system.length)
        function.addParameter('y', dtype='d', direction=function.IN, unit=nbody_system.length)
        function.addParameter('z', dtype='d', direction=function.IN, unit=nbody_system.length)
        function.result_type = 'i'
        return function

    @legacy_function   
    def get_upvector():
        """ specify the orientation of the image by setting the direction vector of image y """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT, unit=units.none)
        function.addParameter('y', dtype='d', direction=function.OUT, unit=units.none)
        function.addParameter('z', dtype='d', direction=function.OUT, unit=units.none)
        function.result_type = 'i'
        return function

    @legacy_function   
    def get_image_angle():
        """ angle of image in x direction (for perpective proj.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_angle', dtype='d', direction=function.OUT, unit=units.deg)
        function.result_type = 'i'
        return function

    @legacy_function   
    def set_image_angle():
        """ angle of image in x direction (for perpective proj.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_angle', dtype='d', direction=function.IN, unit=units.deg)
        function.result_type = 'i'
        return function

    @legacy_function   
    def get_image_ratio():
        """ width/height of image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_ratio', dtype='d', direction=function.OUT, unit=None)
        function.result_type = 'i'
        return function

    @legacy_function
    def viewer():
        function = LegacyFunctionSpecification()  
        return function

    @legacy_function
    def trigger_partremoval():
        function = LegacyFunctionSpecification()  
        return function

    @legacy_function
    def trigger_viewer_refresh():
        function = LegacyFunctionSpecification()  
        return function

    def start_viewer(self):
        self.viewer()
        
class FiDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_doc+"\n\n"+instance.parameters.__doc__

class Fi(GravitationalDynamics, GravityFieldCode):
    
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

    def define_properties(self, handler):
        GravitationalDynamics.define_properties(self, handler)
        handler.add_property("get_thermal_energy")
        handler.add_property("get_total_energy")
    
    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)
        
        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        handler.add_method('END', 'initialize_code')

        handler.add_method('EDIT', 'new_dm_particle')
        handler.add_method('UPDATE', 'new_dm_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_dm_particle', False)
        handler.add_method('EDIT', 'new_sph_particle')
        handler.add_method('UPDATE', 'new_sph_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        handler.add_method('EDIT', 'new_star_particle')
        handler.add_method('UPDATE', 'new_star_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_star_particle', False)
        handler.add_method('RUN', 'get_velocity')
        handler.add_method('RUN', 'get_acceleration')
        handler.add_method('RUN', 'get_internal_energy')
        handler.add_method('RUN', 'get_dinternal_energy_dt')
        handler.add_method('RUN', 'get_smoothing_length')
        handler.add_method('RUN', 'get_density')
        handler.add_method('RUN', 'get_pressure')
        handler.add_method('RUN', 'get_star_tform')
        handler.add_method('RUN', 'get_state_sph')
        handler.add_method('RUN', 'get_state_star')
        
        handler.remove_transition('EVOLVED', 'RUN', 'synchronize_model')
        handler.add_transition('EVOLVED', 'UPDATED', 'update_particle_set')
        handler.add_transition('UPDATED', 'RUN', 'synchronize_model')
        
        handler.add_method('RUN', 'get_kinetic_energy')
        handler.add_method('RUN', 'get_potential_energy')
        handler.add_method('RUN', 'get_thermal_energy')
        handler.add_method('RUN', 'get_total_energy')
        handler.add_method('RUN', 'get_total_radius')
        handler.add_method('RUN', 'get_center_of_mass_position')
        handler.add_method('RUN', 'get_center_of_mass_velocity')
        handler.add_method('RUN', 'get_total_mass')
        handler.add_method('RUN', 'get_time')
        handler.add_method('EDIT', 'get_time')
        handler.add_method('UPDATE', 'get_time')
        handler.add_method('INITIALIZED', 'get_time')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'get_time')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'get_time')
        handler.add_method('CHANGE_PARAMETERS_UPDATE', 'get_time')
        
        
        handler.add_method('RUN', 'get_hydro_state_at_point')

        handler.add_method('EDIT', 'get_gravity_at_point')
        handler.add_method('EDIT', 'get_potential_at_point')
        
        
        self.stopping_conditions.define_state(handler)
    

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_eps2", 
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        
        handler.add_method_parameter(
            "get_dtime", 
            "set_dtime",
            "timestep", 
            "timestep for system", 
            default_value = 1.0 | nbody_system.time
        ) 
        
        
        handler.add_boolean_parameter(
            "get_radiate",
            "set_radiate",
            "radiation_flag",
            "Radiation flag. True means: radiation (i.e. radiative cooling/heating) is included. "
                "False means: no radiation, and implies no star formation.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_starform",
            "set_starform",
            "star_formation_flag",
            "Star-formation flag. True means: star formation is included. "
                "False means: no star formation included.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_use_hydro",
            "set_use_hydro",
            "use_hydro_flag",
            "Hydrodynamics flag. True means: SPH hydro included, False means: gravity only.",
            True
        )
        
        handler.add_boolean_parameter(
            "get_sqrttstp",
            "set_sqrttstp",
            "square_root_timestep_flag",
            "Square-root-timestep flag. True means: use sqrt(eps/acc) timestep criterion.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_acc_tstp",
            "set_acc_tstp",
            "acc_timestep_flag",
            "Acceleration-timestep flag. True means: use vref/acc timestep criterion.",
            True
        )
        
        handler.add_boolean_parameter(
            "get_freetstp",
            "set_freetstp",
            "freeform_timestep_flag",
            "Freeform-timestep flag. True means: use freeform (v/freev)**freevexp * (a/freea)**freeaexp timestep criterion.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_usequad",
            "set_usequad",
            "quadrupole_moments_flag",
            "Quadrupole-moments flag. True means: calculate and use quadrupole cell moments.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_directsum",
            "set_directsum",
            "direct_sum_flag",
            "Direct-summation flag. True means: direct N**2 gravity summation.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_selfgrav",
            "set_selfgrav",
            "self_gravity_flag",
            "Self-gravity flag. False means: self-gravity is not used, only external potentials.",
            True
        )
        
        handler.add_boolean_parameter(
            "get_fixthalo",
            "set_fixthalo",
            "fixed_halo_flag",
            "Fixed-halo flag. True means: use fixed (spherical) potential.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_adaptive_eps",
            "set_adaptive_eps",
            "adaptive_smoothing_flag",
            "Adaptive-smoothing flag. True means: use of adaptive gravity smoothing for all particles.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_gdgop",
            "set_gdgop",
            "gadget_cell_opening_flag",
            "Gadget-cell-opening flag. True means: use of Gadget cell opening criterion.",
            True
        )
        
        handler.add_boolean_parameter(
            "get_smoothinput",
            "set_smoothinput",
            "smooth_input_flag",
            "Smooth-input flag. True means: smooth input SPH properties.",
            False
        )
        
        handler.add_boolean_parameter(
            "get_consph",
            "set_consph",
            "conservative_sph_flag",
            "Conservative-SPH flag. True means: use Springel & Hernquist conservative SPH form (currently the only option).",
            True
        )
        
        handler.add_boolean_parameter(
            "get_sphinit",
            "set_sphinit",
            "sph_dens_init_flag",
            "SPH-density-init flag. True means: initialize sph density and h_smooth (most probably useless for AMUSE interface).",
            True
        )
        
        handler.add_boolean_parameter(
            "get_uentropy",
            "set_uentropy",
            "integrate_entropy_flag",
            "Integrate-entropy flag. True means: integrate entropy, else: internal energy.",
            True
        )
        
        handler.add_boolean_parameter(
            "get_isotherm",
            "set_isotherm",
            "isothermal_flag",
            "Isothermal flag. True means: isothermal gas (requires integrate_entropy_flag == False).",
            False
        )
        
        handler.add_boolean_parameter(
            "get_eps_is_h",
            "set_eps_is_h",
            "eps_is_h_flag",
            "Eps-is-h flag. True means: set gas particles gravitational epsilon to h (SPH smoothing length).",
            True
        )

        handler.add_boolean_parameter(
            "get_balsara",
            "set_balsara",
            "balsara_flag",
            "balsara flag. True means: use Balsara viscosity limiter.",
            False
        )

        handler.add_boolean_parameter(
            "get_mingaseps",
            "set_mingaseps",
            "enforce_min_sph_grav_softening_flag",
            "mingaseps flag. True means: enforce minimum gas grav eps.",
            False
        )
        
        
        handler.add_method_parameter(
            "get_firstsnap", 
            "set_firstsnap",
            "first_snapshot", 
            "The number of the first snapshot.", 
            default_value = 0
        )
        
        handler.add_method_parameter(
            "get_stepout", 
            "set_stepout",
            "output_interval", 
            "The number of steps between output.", 
            default_value = 5
        )
        
        handler.add_method_parameter(
            "get_steplog", 
            "set_steplog",
            "log_interval", 
            "The number of steps between logs.", 
            default_value = 5
        )
        
        handler.add_method_parameter(
            "get_max_tbin", 
            "set_max_tbin",
            "maximum_time_bin", 
            "The maximum time bin (dtime*2**-max_tbin=minimum time step).", 
            default_value = 4096
        )
        
        handler.add_method_parameter(
            "get_minppbin", 
            "set_minppbin",
            "minimum_part_per_bin", 
            "The minimum number of particles per time bin.", 
            default_value = 1
        )
        
        handler.add_method_parameter(
            "get_targetnn", 
            "set_targetnn",
            "targetnn", 
            "The target number of neighbour particles for variable gravitational eps.", 
            default_value = 32
        )
        
        handler.add_method_parameter(
            "get_verbosity", 
            "set_verbosity",
            "verbosity", 
            "The level of terminal output (0=minimum).", 
            default_value = 0
        )
        
        handler.add_method_parameter(
            "get_nsmooth", 
            "set_nsmooth",
            "n_smooth", 
            "The target number of SPH neighbours.", 
            default_value = 64
        )
        
        
        handler.add_method_parameter(
            "get_pboxsize", 
            "set_pboxsize",
            "periodic_box_size", 
            "The size of simulation domain box (particles outside get deleted).", 
            default_value = 10000.0 | nbody_system.length
        )
        
        handler.add_method_parameter(
            "get_unitm_in_msun", 
            "set_unitm_in_msun",
            "code_mass_unit", 
            "The code mass unit (in Msun, 10^9 standard).", 
            default_value = 1.0e9 | units.MSun
        )
        
        handler.add_method_parameter(
            "get_unitl_in_kpc", 
            "set_unitl_in_kpc",
            "code_length_unit", 
            "The code length unit (in kpc, 1 standard).", 
            default_value = 1.0 | units.kpc
        )
        
        handler.add_method_parameter(
            "get_tstepcrit", 
            "set_tstepcrit",
            "sqrt_timestep_crit_constant", 
            "Square-root-timestep criterion constant (unitless,standard=1.).", 
            default_value = 1.0
        )
        
        handler.add_method_parameter(
            "get_tstpcr2", 
            "set_tstpcr2",
            "acc_timestep_crit_constant", 
            "Acceleration-timestep criterion constant (unitless,standard=0.25).", 
            default_value = 0.25
        )
        
        handler.add_method_parameter(
            "get_freev", 
            "set_freev",
            "free_timestep_crit_constant_v", 
            "Freeform-timestep criterion constant v.", 
            default_value = 0.5
        )
        
        handler.add_method_parameter(
            "get_freea", 
            "set_freea",
            "free_timestep_crit_constant_a", 
            "Freeform-timestep criterion constant a.", 
            default_value = 0.35
        )
        
        handler.add_method_parameter(
            "get_freevexp", 
            "set_freevexp",
            "free_timestep_crit_constant_vexp", 
            "Freeform-timestep criterion constant v_exp.", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_freeaexp", 
            "set_freeaexp",
            "free_timestep_crit_constant_aexp", 
            "Freeform-timestep criterion constant a_exp.", 
            default_value = -1.0
        )
        
        handler.add_method_parameter(
            "get_bh_tol", 
            "set_bh_tol",
            "opening_angle", 
            "Opening angle, theta, for building the tree: between 0 and 1 (unitless, 0.5).", 
            default_value = 0.5
        )
        
        handler.add_method_parameter(
            "get_gdgtol", 
            "set_gdgtol",
            "gadget_cell_opening_constant", 
            "Gadget-cell-openings criterion parameter  (unitless, .01)", 
            default_value = 0.01
        )
        
        handler.add_method_parameter(
            "get_nn_tol", 
            "set_nn_tol",
            "nn_tol", 
            "The fractional tolerance in nn_target (0.1).", 
            default_value = 0.1
        )
        
        handler.add_method_parameter(
            "get_epsgas", 
            "set_epsgas",
            "gas_epsilon", 
            "The gas gravitational smoothing epsilon.", 
            default_value = 0.005 | nbody_system.length
        )
        
        handler.add_method_parameter(
            "get_gamma", 
            "set_gamma",
            "gamma", 
            "gas polytropic index (1.6666667)", 
            default_value = 1.6666667
        )
        
        handler.add_method_parameter(
            "get_alpha", 
            "set_alpha",
            "artificial_viscosity_alpha", 
            "SPH artificial viscosity alpha parameter (0.5)", 
            default_value = 0.5
        )
        
        handler.add_method_parameter(
            "get_beta", 
            "set_beta",
            "beta", 
            "SPH artificial viscosity beta parameter (2*alpha=1.0)", 
            default_value = 1.0
        )
        
        handler.add_method_parameter(
            "get_epssph", 
            "set_epssph",
            "sph_artificial_viscosity_eps", 
            "SPH artificial viscosity safety against divergence (0.01)", 
            default_value = 0.01
        )
        
        handler.add_method_parameter(
            "get_courant", 
            "set_courant",
            "courant", 
            "SPH courant condition parameter (0.3)", 
            default_value = 0.3
        )
        
        handler.add_method_parameter(
            "get_removgas", 
            "set_removgas",
            "min_gas_part_mass", 
            "minimum gas particle mass (fraction of initial (average) mass)", 
            default_value = 0.25
        )
        
        handler.add_method_parameter(
            "get_consthsm", 
            "set_consthsm",
            "sph_h_const", 
            "SPH smoothing length if constant", 
            default_value = 0.2 | nbody_system.length
        )
        
        handler.add_method_parameter(
            "get_nsmtol", 
            "set_nsmtol",
            "n_smooth_tol", 
            "fractional tolerance in number of SPH neighbours", 
            default_value = 0.1
        )
        
        handler.add_method_parameter(
            "get_graineff", 
            "set_graineff",
            "grain_heat_eff", 
            "FUV grain heating efficiency parameter (unitless, 0.05)", 
            default_value = 0.05
        )
        
        handler.add_method_parameter(
            "get_crionrate", 
            "set_crionrate",
            "zeta_cr_ion_rate", 
            "primary cosmic ray ionization rate, zeta (in units of 1.8e-17 sec^-1, 1.)", 
            default_value = 3.6 | 1.8e-17 * units.s**-1
        )
        
        handler.add_method_parameter(
            "get_heat_par1", 
            "set_heat_par1",
            "heat_par1", 
            "additional heating 1 parameter (0.0)", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_heat_par2", 
            "set_heat_par2",
            "heat_par2", 
            "additional heating 2 parameter (0.0)", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_cool_par", 
            "set_cool_par",
            "cool_par", 
            "additional cooling parameter (1.0)", 
            default_value = 1.0
        )
        
        handler.add_method_parameter(
            "get_optdepth", 
            "set_optdepth",
            "optical_depth", 
            "1/(mean free path) for UV photons (code length **-1, 0.0)", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_tcollfac", 
            "set_tcollfac",
            "star_form_delay_fac", 
            "star formation delay parameter (unitless, 1)", 
            default_value = 1.0
        )
        
        handler.add_method_parameter(
            "get_masscrit", 
            "set_masscrit",
            "star_form_mass_crit", 
            "star formation cloud reference mass (Msun, 1.e5)", 
            default_value = 1.0e5 | units.MSun
        )
        
        handler.add_method_parameter(
            "get_sfeff", 
            "set_sfeff",
            "star_form_eff", 
            "gas particle mass fraction converted to stars (0.125)", 
            default_value = 0.25
        )
        
        handler.add_method_parameter(
            "get_tbubble", 
            "set_tbubble",
            "supernova_duration", 
            "Supernova activity time, (Myr, 3.e7)", 
            default_value = 3.0e7 | units.Myr
        )
        
        handler.add_method_parameter(
            "get_sne_eff", 
            "set_sne_eff",
            "supernova_eff", 
            "Supernova feedback coupling efficiency, (0.0)", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_tsnbeg", 
            "set_tsnbeg",
            "t_supernova_start", 
            "Supernova feedback start time, (Myr, 3.e6)", 
            default_value = 3.0e6 | units.Myr
        )
        
        handler.add_method_parameter(
            "get_rhomax", 
            "set_rhomax",
            "max_density", 
            "Maximum density in case of star formation (force SF if exceeded, ignored if star formation is off)", 
            default_value = 100.0 |nbody_system.density
        )
        
        handler.add_method_parameter(
            "get_halofile", 
            "set_halofile",
            "halofile", 
            "Path to initial halo model file, relative to the Fi data directory (none)", 
            default_value = "none"
        )
        
        handler.add_method_parameter(
            "get_feedback", 
            "set_feedback",
            "feedback", 
            "feedback model (fuv, pres, kine, solo, solh)", 
            default_value = "fuv"
        )
        
        handler.add_method_parameter(
            "get_sfmode", 
            "set_sfmode",
            "star_formation_mode", 
            "star formation model (gerritsen, nieuw)", 
            default_value = "gerritsen"
        )
        
        handler.add_method_parameter(
            "get_hupdatemethod", 
            "set_hupdatemethod",
            "h_update_method", 
            "SPH smoothing length criterion (at the moment always 'mass')", 
            default_value = "mass"
        )
        
        handler.add_method_parameter(
            "get_sph_visc", 
            "set_sph_visc",
            "sph_viscosity", 
            "SPH viscosity (sph,sphv, bulk). Note: not all may work.", 
            default_value = "sph"
        )
        
        handler.add_method_parameter(
            "get_fi_data_directory", 
            "set_fi_data_directory",
            "fi_data_directory", 
            "Name of the Fi data directory", 
            default_value = ""
        )

        handler.add_boolean_parameter(
            "get_periodic_boundaries_flag",
            None,
            "periodic_boundaries_flag",
            "Periodic boundaries flag. True means: use periodic boundary conditions (read-only)",
            False
        )
        
        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )
        
        self.stopping_conditions.define_parameters(handler)        
    
    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['dm_particles','gas_particles','star_particles'], 
            index_to_default_set = 0)
        
        handler.define_set('dm_particles', 'id')
        handler.set_new('dm_particles', 'new_dm_particle')
        handler.set_delete('dm_particles', 'delete_particle')
        handler.add_setter('dm_particles', 'set_state')
        handler.add_getter('dm_particles', 'get_state')
        handler.add_setter('dm_particles', 'set_mass')
        handler.add_getter('dm_particles', 'get_mass', names = ('mass',))
        handler.add_setter('dm_particles', 'set_position')
        handler.add_getter('dm_particles', 'get_position')
        handler.add_setter('dm_particles', 'set_radius')
        handler.add_getter('dm_particles', 'get_radius')
        handler.add_setter('dm_particles', 'set_velocity')
        handler.add_getter('dm_particles', 'get_velocity')
        
        handler.define_set('gas_particles', 'id')
        handler.set_new('gas_particles', 'new_sph_particle')
        handler.set_delete('gas_particles', 'delete_particle')
        handler.add_setter('gas_particles', 'set_state_sph')
        handler.add_getter('gas_particles', 'get_state_sph')
        handler.add_setter('gas_particles', 'set_mass')
        handler.add_getter('gas_particles', 'get_mass', names = ('mass',))
        handler.add_getter('gas_particles', 'get_radius')
        handler.add_setter('gas_particles', 'set_radius')
        handler.add_setter('gas_particles', 'set_position')
        handler.add_getter('gas_particles', 'get_position')
        handler.add_setter('gas_particles', 'set_velocity')
        handler.add_getter('gas_particles', 'get_velocity')
        handler.add_setter('gas_particles', 'set_internal_energy')
        handler.add_getter('gas_particles', 'get_internal_energy')
        handler.add_getter('gas_particles', 'get_dinternal_energy_dt')
        handler.add_setter('gas_particles', 'set_smoothing_length')
        handler.add_getter('gas_particles', 'get_smoothing_length')
        handler.add_getter('gas_particles', 'get_density', names = ('rho',))
        handler.add_getter('gas_particles', 'get_density', names = ('density',))
        handler.add_getter('gas_particles', 'get_pressure')
        
        handler.define_set('star_particles', 'id')
        handler.set_new('star_particles', 'new_star_particle')
        handler.set_delete('star_particles', 'delete_particle')
        handler.add_setter('star_particles', 'set_state_star')
        handler.add_getter('star_particles', 'get_state_star')
        handler.add_setter('star_particles', 'set_mass')
        handler.add_getter('star_particles', 'get_mass', names = ('mass',))
        handler.add_setter('star_particles', 'set_position')
        handler.add_getter('star_particles', 'get_position')
        handler.add_setter('star_particles', 'set_radius')
        handler.add_getter('star_particles', 'get_radius')
        handler.add_setter('star_particles', 'set_velocity')
        handler.add_getter('star_particles', 'get_velocity')
        handler.add_setter('star_particles', 'set_star_tform')
        handler.add_getter('star_particles', 'get_star_tform')
        
        self.stopping_conditions.define_particle_set(handler, 'particles')

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        handler.add_method(
            "set_velocity",
            (
                handler.INDEX,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_velocity",
            (
                handler.INDEX,
            ),
            (
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
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
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
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
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state_sph",
            (
                handler.INDEX,
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
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_state_sph",
            (
                handler.INDEX,
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
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_internal_energy",
            (
                handler.INDEX,
                nbody_system.specific_energy,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_internal_energy",
            (
                handler.INDEX,
            ),
            (
                nbody_system.specific_energy,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_dinternal_energy_dt",
            (
                handler.INDEX,
            ),
            (
                nbody_system.specific_energy/nbody_system.time,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_smoothing_length",
            (handler.INDEX, nbody_system.length),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_smoothing_length",
            (handler.INDEX,),
            (nbody_system.length, handler.ERROR_CODE)
        )
        handler.add_method(
            "get_density",
            (handler.INDEX,),
            (nbody_system.density, handler.ERROR_CODE)
        )
        handler.add_method(
            "get_pressure",
            (handler.INDEX,),
            (nbody_system.pressure, handler.ERROR_CODE)
        )
        
        handler.add_method(
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
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state_star",
            (
                handler.INDEX,
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
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_state_star",
            (
                handler.INDEX,
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
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "set_star_tform",
            (
                handler.INDEX,
                nbody_system.time,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_star_tform",
            (
                handler.INDEX,
            ),
            (
                nbody_system.time,
                handler.ERROR_CODE
            )
        )
        
        handler.add_method(
            'get_hydro_state_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length,
                nbody_system.speed, nbody_system.speed, nbody_system.speed),
            (nbody_system.density, nbody_system.momentum_density, nbody_system.momentum_density, 
                nbody_system.momentum_density, nbody_system.energy_density, handler.ERROR_CODE)
        )
        
        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_dtime",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_dtime",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_firstsnap",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_firstsnap",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_stepout",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_stepout",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_steplog",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_steplog",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_max_tbin",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_max_tbin",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_minppbin",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_minppbin",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_targetnn",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_targetnn",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_verbosity",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_verbosity",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_nsmooth",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_nsmooth",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_pboxsize",
            (),
            (nbody_system.length, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_pboxsize",
            (nbody_system.length, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_unitm_in_msun",
            (),
            (units.MSun, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_unitm_in_msun",
            (units.MSun, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_unitl_in_kpc",
            (),
            (units.kpc, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_unitl_in_kpc",
            (units.kpc, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_tstepcrit",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_tstepcrit",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_tstpcr2",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_tstpcr2",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_freev",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_freev",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_freea",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_freea",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_freevexp",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_freevexp",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_freeaexp",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_freeaexp",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_bh_tol",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_bh_tol",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_gdgtol",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_gdgtol",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_nn_tol",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_nn_tol",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_epsgas",
            (),
            (nbody_system.length, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_epsgas",
            (nbody_system.length, ),
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
            "get_alpha",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_alpha",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_beta",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_beta",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_epssph",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_epssph",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_courant",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_courant",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_removgas",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_removgas",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_consthsm",
            (),
            (nbody_system.length, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_consthsm",
            (nbody_system.length, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_nsmtol",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_nsmtol",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_graineff",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_graineff",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_crionrate",
            (),
            (1.8e-17 * units.s**-1, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_crionrate",
            (1.8e-17 * units.s**-1, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_heat_par1",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_heat_par1",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_heat_par2",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_heat_par2",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_cool_par",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_cool_par",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_optdepth",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_optdepth",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_tcollfac",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_tcollfac",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_masscrit",
            (),
            (units.MSun, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_masscrit",
            (units.MSun, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_sfeff",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_sfeff",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_tbubble",
            (),
            (units.Myr, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_tbubble",
            (units.Myr, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_sne_eff",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_sne_eff",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_tsnbeg",
            (),
            (units.Myr, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_tsnbeg",
            (units.Myr, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_rhomax",
            (),
            (nbody_system.density, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_rhomax",
            (nbody_system.density, ),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_halofile",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_halofile",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_feedback",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_feedback",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_sfmode",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_sfmode",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_hupdatemethod",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_hupdatemethod",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_sph_visc",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_sph_visc",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_fi_data_directory",
            (),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_fi_data_directory",
            (handler.NO_UNIT,),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_thermal_energy",
            (),
            (nbody_system.energy, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_total_energy",
            (),
            (nbody_system.energy, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_number_of_sph_particles_removed",
            (
            ),
            (
                handler.NO_UNIT,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_id_of_removed_sph_particle",
            (
                handler.NO_UNIT,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE
            )
        )

        self.stopping_conditions.define_methods(handler)       
        
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
        
        indices_in_update_list = list(range(number_of_removed_particles))
        indices_to_remove = self.get_id_of_removed_sph_particle(indices_in_update_list)
        
        incode_storage = self.gas_particles._remove_indices_in_attribute_storage(indices_to_remove)

    def get_total_energy(self):
        return self.__getattr__('get_total_energy')()
        
class FiViewer(Fi):

    def __init__(self, convert_nbody = None, mode = 'normal', **options):
        Fi.__init__(self, convert_nbody = convert_nbody, mode = mode, use_gl = True , **options)

    def define_parameters(self, handler):
        Fi.define_parameters(self,handler)

        handler.add_method_parameter(
            "get_viewpoint",
            "set_viewpoint",
            "viewpoint",
            "viewpoint (location of the camera)",
            [0,1,0] | nbody_system.length, is_vector=True
        )

        handler.add_method_parameter(
            "get_image_target",
            "set_image_target",
            "image_target",
            "image_target (location the camera points to)",
            [0,0,0] | nbody_system.length, is_vector=True
        )
        
        handler.add_method_parameter(
            "get_upvector",
            None,
            "upvector",
            "upvector (which direction points up, readonly)",
            [0,0,1] , is_vector=True
        )        

        handler.add_method_parameter(
            "get_image_angle", 
            "set_image_angle",
            "image_angle", 
            "image angle - vertical (?)!! ", 
            default_value = 45 | units.deg
        )
        
        handler.add_method_parameter(
            "get_image_ratio", 
            None,
            "image_ratio", 
            "image width/height (readonly)", 
            default_value = 1.
        )        
        
        
    def define_state(self, handler):
        handler.set_initial_state('UNINITIALIZED')
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        handler.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        handler.add_transition('END', 'STOPPED', 'stop', False)
        handler.add_method('STOPPED', 'stop')
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_method('EDIT', 'new_dm_particle')
        handler.add_method('EDIT', 'new_sph_particle')
        handler.add_method('EDIT', 'new_star_particle')
        handler.add_method('UPDATE', 'delete_particle')
        handler.add_transition('EDIT', 'UPDATE', 'delete_particle')
        handler.add_transition('UPDATE', 'EDIT', 'trigger_partremoval')
                  
class FiMapInterface(CodeInterface):   

    use_modules=['StoppingConditions','MapInterface']
    
    MODE_NORMAL = 'normal'
    MODE_NORMAL_OPENMP = 'openmp'
    
    def __init__(self, mode = MODE_NORMAL,  **options):
        self.mode = mode
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(mode), **options)
                     
    def name_of_the_worker(self, mode):
        if mode == self.MODE_NORMAL:
            return 'fi_worker_map'
        elif mode == self.MODE_NORMAL_OPENMP:
            return 'fi_worker_map_mp'
        else:
            return 'fi_worker_map'
  
    @legacy_function    
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def cleanup_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

  
    @legacy_function    
    def generate_projection():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def init_map():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def reset_map():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def erase_map():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('weight', dtype='d', direction=function.IN, default = 1.)
        function.addParameter('radius', dtype='d', direction=function.IN, default = 0.)
        function.addParameter('opacity_area', dtype='d', direction=function.IN, default = 0.)
        function.addParameter('index', dtype='i', direction=function.IN, default = -1)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('weight', dtype='d', direction=function.IN, default = 1.)
        function.addParameter('radius', dtype='d', direction=function.IN, default = 0.)
        function.addParameter('opacity_area', dtype='d', direction=function.IN, default = 0.)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_weight():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('weight', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_opacity_area():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('opacity_area', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('weight', dtype='d', direction=function.OUT)
        function.addParameter('radius', dtype='d', direction=function.OUT)
        function.addParameter('opacity_area', dtype='d', direction=function.OUT)
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
    def set_random_seed():
        """ random seed to use """            
        function = LegacyFunctionSpecification()  
        function.addParameter('random_seed', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_random_seed():
        """ random seed to use """            
        function = LegacyFunctionSpecification()  
        function.addParameter('random_seed', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_extinction_flag():
        """ whether to use the particle opacities (0=no, 1= yes) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('extinction_flag', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_extinction_flag():
        """ whether to use the particle opacities (0=no, 1= yes) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('extinction_flag', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function
    def set_minimum_distance():
        """ minimum distance to the camera particles can have """            
        function = LegacyFunctionSpecification()  
        function.addParameter('minimum_distance', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_minimum_distance():
        """ minimum distance to the camera particles can have """            
        function = LegacyFunctionSpecification()  
        function.addParameter('minimum_distance', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_image_angle():
        """ angle of image in x direction (for perpective proj.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_angle', dtype='d', direction=function.IN,unit=units.deg)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_image_angle():
        """ angle of image in x direction (for perpective proj.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_angle', dtype='d', direction=function.OUT,unit=units.deg)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_image_width():
        """ width of image in x direction (for parallel proj.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_width', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_image_width():
        """ angle of image in x direction (for parallel proj.) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('image_width', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_image_pixel_size():
        """ pixel size of generated image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('nx', dtype='i', direction=function.IN)
        function.addParameter('ny', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_image_pixel_size():
        """ pixel size of generated image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('nx', dtype='i', direction=function.OUT)
        function.addParameter('ny', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    def get_index_range_inclusive(self):
        ni,nj,error = self.get_image_pixel_size()
        return (1, ni, 1, nj)


    @legacy_function
    def set_image_target():
        """ target point of image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.IN)
        function.addParameter('y', dtype='d', direction=function.IN)
        function.addParameter('z', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_image_target():
        """ target point of image """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT)
        function.addParameter('y', dtype='d', direction=function.OUT)
        function.addParameter('z', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_viewpoint():
        """ camera position (for perspective proj) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.IN)
        function.addParameter('y', dtype='d', direction=function.IN)
        function.addParameter('z', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_viewpoint():
        """ camera position (for perspective proj) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT)
        function.addParameter('y', dtype='d', direction=function.OUT)
        function.addParameter('z', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_projection_direction():
        """ direction of projection (for parallel proj) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.IN)
        function.addParameter('y', dtype='d', direction=function.IN)
        function.addParameter('z', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_projection_direction():
        """ direction of projection (for parallel proj) """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT)
        function.addParameter('y', dtype='d', direction=function.OUT)
        function.addParameter('z', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_upvector():
        """ specify the orientation of the image by setting the direction vector of image y """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.IN)
        function.addParameter('y', dtype='d', direction=function.IN)
        function.addParameter('z', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_upvector():
        """ specify the orientation of the image by setting the direction vector of image y """            
        function = LegacyFunctionSpecification()  
        function.addParameter('x', dtype='d', direction=function.OUT)
        function.addParameter('y', dtype='d', direction=function.OUT)
        function.addParameter('z', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_projection_mode():
        """ projection mode (parallel or projection """
        function = LegacyFunctionSpecification()  
        function.addParameter('projection_mode', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function   
    def get_projection_mode():
        """ projection mode (parallel or projection """
        function = LegacyFunctionSpecification()  
        function.addParameter('projection_mode', dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function


    @legacy_function    
    def get_image():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['pixel_value']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_opdepth_map():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['pixel_opacity']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

class FiMap(CommonCode):

    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        
        CommonCode.__init__(self,  FiMapInterface(**options), **options)
    
    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())

    def define_methods(self, handler):
        handler.add_method(
             'new_particle', 
              (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                handler.NO_UNIT,
                generic_unit_system.length,
                generic_unit_system.length**2,
                handler.NO_UNIT,                
              ), 
              (
                handler.INDEX,                    
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'set_state', 
              (
                handler.INDEX,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                handler.NO_UNIT,
                generic_unit_system.length,
                generic_unit_system.length**2,
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'set_weight', 
              (
                handler.INDEX,                    
                handler.NO_UNIT,
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'set_opacity_area', 
              (
                handler.INDEX,                    
                generic_unit_system.length**2,
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'get_state', 
              (
                handler.INDEX,
              ), 
              (
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                handler.NO_UNIT,
                generic_unit_system.length,
                generic_unit_system.length**2,
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'delete_particle', 
              (
                handler.INDEX,                    
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        
        handler.add_method(
             'set_minimum_distance', 
              (
                generic_unit_system.length,                    
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'get_minimum_distance', 
              (
              ), 
              (
                generic_unit_system.length,                    
                handler.ERROR_CODE,
              )
        )

        handler.add_method(
             'set_image_width', 
              (
                generic_unit_system.length,                    
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'get_image_width', 
              (
              ), 
              (
                generic_unit_system.length,                    
                handler.ERROR_CODE,
              )
        )

        handler.add_method(
             'set_image_target', 
              (
                generic_unit_system.length,                    
                generic_unit_system.length,                    
                generic_unit_system.length,                    
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'get_image_target', 
              (
              ), 
              (
                generic_unit_system.length,                    
                generic_unit_system.length,                    
                generic_unit_system.length,                    
                handler.ERROR_CODE,
              )
        )

        handler.add_method(
             'set_viewpoint', 
              (
                generic_unit_system.length,                    
                generic_unit_system.length,                    
                generic_unit_system.length,                    
              ), 
              (
                handler.ERROR_CODE,
              )
        )
        handler.add_method(
             'get_viewpoint', 
              (
              ), 
              (
                generic_unit_system.length,                    
                generic_unit_system.length,                    
                generic_unit_system.length,                    
                handler.ERROR_CODE,
              )
        )

        handler.add_method(
             'get_image', 
              (
                handler.INDEX,
                handler.INDEX,
              ), 
              (
                handler.NO_UNIT,
                handler.ERROR_CODE,
              )
        )

        handler.add_method(
             'get_opdepth_map', 
              (
                handler.INDEX,
                handler.INDEX,
              ), 
              (
                handler.NO_UNIT,
                handler.ERROR_CODE,
              )
        )


    def define_particle_sets(self, handler):
        handler.define_grid('image')
        handler.set_grid_range('image', 'get_index_range_inclusive')    
        handler.add_getter('image', 'get_image')
        handler.add_getter('image', 'get_opdepth_map')

        handler.define_set('particles', 'id')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_state')
        handler.add_setter('particles', 'set_weight')
        handler.add_setter('particles', 'set_opacity_area')
        handler.add_getter('particles', 'get_state')


    def define_parameters(self, handler):
        
        handler.add_method_parameter(
            "get_projection_mode", 
            "set_projection_mode",
            "projection_mode", 
            "projection mode (parallel or perspective)", 
            default_value = "parallel"
        )


        handler.add_method_parameter(
            "get_random_seed", 
            "set_random_seed",
            "random_seed", 
            "random seed (used for clusters)", 
            default_value = 45678910
        )

        handler.add_method_parameter(
            "get_minimum_distance", 
            "set_minimum_distance",
            "minimum_distance", 
            "near clipping plane", 
            default_value = 0.001 | generic_unit_system.length
        )

        handler.add_boolean_parameter(
            "get_extinction_flag", 
            "set_extinction_flag",
            "extinction_flag", 
            "extinction flag (whether to use extinction)", 
            default_value = False
        )

        handler.add_method_parameter(
            "get_image_angle", 
            "set_image_angle",
            "image_angle", 
            "image angle - horizontal (used in perspective proj.)", 
            default_value = 45 | units.deg
        )

        handler.add_method_parameter(
            "get_image_width", 
            "set_image_width",
            "image_width", 
            "image width - horizontal (used in parallel proj.)", 
            default_value = 1 | generic_unit_system.length
        )

        
        handler.add_caching_parameter(
            "set_image_pixel_size", 
            "nx",
            "nx", 
            "image pixel size (horizontal)",
            640,
        )
        handler.add_caching_parameter(
            "set_image_pixel_size", 
            "ny",
            "ny", 
            "image pixel size (vertical)",
            480,
        )
        handler.add_vector_parameter(
            "image_size",
            "image pixel size",
            ("nx", "ny")
        )

        handler.add_caching_parameter(
            "set_image_target", 
            "x",
            "target_x", 
            "x coordinate of the point which the image centers on",
            0 | generic_unit_system.length,
        )
        handler.add_caching_parameter(
            "set_image_target", 
            "y",
            "target_y", 
            "y coordinate of the point which the image centers on",
            0 | generic_unit_system.length,
        )
        handler.add_caching_parameter(
            "set_image_target", 
            "z",
            "target_z", 
            "z coordinate of the point which the image centers on",
            0 | generic_unit_system.length,
        )
        handler.add_vector_parameter(
            "image_target",
            "point which the image centers on",
            ("target_x", "target_y","target_z")
        )

        handler.add_caching_parameter(
            "set_viewpoint", 
            "x",
            "viewpoint_x", 
            "x coordinate of the view point (camera location)",
            0 | generic_unit_system.length,
        )
        handler.add_caching_parameter(
            "set_viewpoint", 
            "y",
            "viewpoint_y", 
            "y coordinate of the view point (camera location)",
            1. | generic_unit_system.length,
        )
        handler.add_caching_parameter(
            "set_viewpoint", 
            "z",
            "viewpoint_z", 
            "z coordinate of the view point (camera location)",
            0 | generic_unit_system.length,
        )
        handler.add_vector_parameter(
            "viewpoint",
            "viewpoint (location of the camera)",
            ("viewpoint_x", "viewpoint_y","viewpoint_z")
        )

        handler.add_caching_parameter(
            "set_upvector", 
            "x",
            "upvector_x", 
            "x component of the up-direction of the image",
            0,
        )
        handler.add_caching_parameter(
            "set_upvector", 
            "y",
            "upvector_y", 
            "y component of the up-direction of the image",
            0,
        )
        handler.add_caching_parameter(
            "set_upvector", 
            "z",
            "upvector_z", 
            "z component of the up-direction of the image",
            1,
        )
        handler.add_vector_parameter(
            "upvector",
            "direction of the up-vector",
            ("upvector_x", "upvector_y","upvector_z")
        )

        handler.add_caching_parameter(
            "set_projection_direction", 
            "x",
            "projection_direction_x", 
            "x component of projection direction (for parallel projections)",
            0,
        )
        handler.add_caching_parameter(
            "set_projection_direction", 
            "y",
            "projection_direction_y", 
            "y component of projection direction (for parallel projections)",
            -1,
        )
        handler.add_caching_parameter(
            "set_projection_direction", 
            "z",
            "projection_direction_z", 
            "z component of projection direction (for parallel projections)",
            0,
        )
        handler.add_vector_parameter(
            "projection_direction",
            "direction of projection (for parallel projection)",
            ("projection_direction_x", "projection_direction_y","projection_direction_z")
        )

    def get_index_range_inclusive(self):
        """
        Returns the min and max values of indices in each
        direction. The range is inclusive, the min index
        and max index both exist and can be queried.
        The total number of cells in one direction
        is max - min + 1.
        """
        nx, ny = self.get_image_pixel_size()
        return (1, nx, 1, ny)
    
    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
              
    def define_state(self, handler): 
        CommonCode.define_state(self, handler)   
        #handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        
        handler.add_transition('INITIALIZED','PREPROJ','commit_parameters')
        handler.add_transition('PREPROJ','PROJ','init_map')
        handler.add_transition('PROJ','ERASE','erase_map')
        handler.add_transition('ERASE','IMAGE','generate_projection')
        handler.add_transition('IMAGE','PROJ','new_particle',False)        
        handler.add_transition('IMAGE','PROJ','delete_particle',False)        
        handler.add_transition('IMAGE','PROJ','set_state',False)
        handler.add_transition('IMAGE','PROJ','set_weight',False)
        handler.add_transition('IMAGE','PROJ','set_opacity_area',False)
        handler.add_transition('PROJ','INITIALIZED','reset_map')
        handler.add_transition('IMAGE','INITIALIZED','reset_map')
        handler.add_method('IMAGE','get_image')
        handler.add_method('IMAGE','get_opdepth_map')
        handler.add_method('PROJ', 'new_particle')
        handler.add_method('PROJ', 'delete_particle')
        handler.add_method('PROJ', 'set_state')
        handler.add_method('PROJ', 'set_weight')
        handler.add_method('PROJ', 'set_opacity_area')

        handler.add_method('INITIALIZED', 'before_set_parameter')  
        handler.add_method('PROJ', 'before_get_parameter')
        handler.add_method('PROJ', 'get_image_pixel_size')
        handler.add_method('IMAGE', 'get_image_pixel_size')
