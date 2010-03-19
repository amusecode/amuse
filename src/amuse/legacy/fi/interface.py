import numpy
from amuse.legacy.interface.gd import GravitationalDynamicsInterface, GravitationalDynamics
#from amuse.support.units import nbody_system
#from amuse.support.units import units
from amuse.legacy import *

class fi(LegacyInterface,GravitationalDynamicsInterface):   
    get_density_at_point=None
    get_temperature_at_point=None
    get_internalenergy_at_point=None
    get_pressure_at_point=None

    get_potential_energy=None
    get_kinetic_energy=None    
    get_thermal_energy=None
    get_total_radius=None
    get_total_mass=None
    get_center_of_mass_position=None
    get_center_of_mass_velocity=None

    get_indices_of_colliding_particles=None

    get_potential=None
    get_mass=None
    set_mass=None
    set_eps2=None
    get_eps2=None
    set_acceleration=None
    get_acceleration=None
    get_radius=None
    set_radius=None
    get_velocity=None
    set_velocity=None
    set_position=None
    get_position=None

        
    def __init__(self, convert_nbody = None):
        LegacyInterface.__init__(self, name_of_the_worker = 'worker')
                     
    def setup_module(self):
        return self.initialize_code()
        
    def cleanup_module(self):
        return self.cleanup_code()
       
    def initialize_particles(self, ignore):
        return self.commit_particles()
        
    def reinitialize_particles(self):
        return self.recommit_particles()
                
    @legacy_function    
    def new_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_sph_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def new_star_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz','tform']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
             
    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        return function

    @legacy_function    
    def get_state_sph():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        return function

    @legacy_function    
    def get_state_star():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz','tform']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        return function

    @legacy_function    
    def set_state():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        return function

    @legacy_function    
    def set_state_star():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        return function


    @legacy_function    
    def set_state_sph():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
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
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function    

    @legacy_function    
    def evolve():
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
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_potential_at_point():
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('phi', dtype='d', direction=function.OUT)
        function.result_type = 'i'
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
    def set_usesph():
        """ set_usesph([0,1]): SPH hydro if 0, grav only if 1 """
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_usesph():
        """ get_usesph(): SPH hydro if 0, grav only if 1 """
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_radiate():
        """ set_radiate([0,1]): rad cooling/heating if 0, not if 1
             radiate false implies starform false """    
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_radiate():
        """ get_radiate(): rad cooling/heating if 0, not if 1
             radiate false implies starform false """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_starform():
        """ set_starform([0,1]): star formation if 0, not if 1 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_starform():
        """ get_starform(): star formation if 0, not if 1 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
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
        """ set_sqrttstp([0,1]): use sqrt(eps/acc) timestep crit if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sqrttstp():
        """ get_sqrttstp(): use sqrt(eps/acc) timestep crit if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_acc_tstp():
        """ set_acc_tstp([0,1]): use vref/acc timestep crit if 0"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_acc_tstp():
        """ get_acc_tstp(): use vref/acc timestep crit if 0"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_freetstp():
        """ set_freetstp([0,1]): use freeform (v/freev)**freevexp * (a/freea)**freeaexp timestep crit if 0"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_freetstp():
        """ get_freetstp(): use freeform timestep crit if 0"""            
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_usequad():
        """ set_usequad([0,1]): calc. and use quadrupole cell moments if 0"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_usequad():
        """ get_usequad(): calc. and use quadrupole cell moments if 0"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_directsum():
        """ set_directsum([0,1]): direct N**2 grav sum if 0"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_directsum():
        """ get_directsum(): direct N**2 grav sum if 0"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_selfgrav():
        """ set_selfgrav([0,1]): calculate self-gravity if 0
          if set to 1, self gravity is not used, only external potentials"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_selfgrav():
        """ get_selfgrav(): calculate self-gravity if 0 """
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_fixthalo():
        """ set_fixthalo([0,1]): use fixed (spherical) potential if 0 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_fixthalo():
        """ get_fixthalo(): use fixed (spherical) potential if 0 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_adaptive_eps():
        """ set_adaptive_eps([0,1]): use of adaptive grav smoothing for all part if 0 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_adaptive_eps():
        """ get_adaptive_eps(): use of adaptive grav smoothing for all part if 0 """    
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_gdgop():
        """ set_gdgop([0,1]): use of gadget cell opening criterion if 0 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_gdgop():
        """ get_gdgop(): use of gadget cell opening criterion if 0 """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_smoothinput():
        """ set_smoothinput([0,1]): smooth input SPH prop. if 0 
         (not working) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_smoothinput():
        """ get_smoothinput(): smooth input SPH prop. if 0 
         (not working) """        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_consph():
        """ set_consph([0,1]): use springel&Hernquist conservative SPH form. if 0 
          at the moment this is only option"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_consph():
        """ get_consph(): use springel&Hernquist conservative SPH form. if 0 
          at the moment this is only option"""                
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_sphinit():
        """ set_sphinit([0,1]): initialize sph dens and hsmooth if 0 
         most probably useless for AMUSE interface"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_sphinit():
        """ set_sphinit([0,1]): initialize sph dens and hsmooth if 0 
         most probably useless for AMUSE interface"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_uentropy():
        """ set_uentropy([0,1]): integrate entropy if 0, internal energy if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_uentropy():
        """ get_uentropy(): integrate entropy if 0, internal energy if 1"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_isotherm():
        """ set_isotherm([0,1]): isothermal gas if 0
          note that isotherm needs set_uentropy(1) (false)"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_isotherm():
        """ get_isotherm(): isothermal gas if 0
          note that isotherm needs set_uentropy(1) (false)"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function   
    def set_eps_is_h():
        """ set_eps_is_h([0,1]): gas particles grav. eps to SPH h if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_eps_is_h():
        """ get_eps_is_h(): gas particles grav. eps to SPH h if 0"""        
        function = LegacyFunctionSpecification()  
        function.addParameter('zeroiftrue', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;


# integers
    @legacy_function
    def set_firstsnap():
        """ no. of first snapshot """
        function = LegacyFunctionSpecification()  
        function.addParameter('firstsnap', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_firstsnap():
        """ no. of first snapshot """
        function = LegacyFunctionSpecification()  
        function.addParameter('firstsnap', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_stepout():
        """ no. of steps between output """
        function = LegacyFunctionSpecification()  
        function.addParameter('stepout', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_stepout():
        """ no. of steps between output """
        function = LegacyFunctionSpecification()  
        function.addParameter('stepout', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_steplog():
        """ no. of steps between logs """
        function = LegacyFunctionSpecification()  
        function.addParameter('steplog', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_steplog():
        """ no. of steps between logs """
        function = LegacyFunctionSpecification()  
        function.addParameter('steplog', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_max_tbin():
        """ maximum time bin (dtime*2**-max_tbin=minimum time step)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('max_tbin', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_max_tbin():
        """ maximum time bin (dtime*2**-max_tbin=minimum time step)"""
        function = LegacyFunctionSpecification()  
        function.addParameter('max_tbin', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function;

    @legacy_function
    def set_minppbin():
        """ target no. of particles per time bin"""
        function = LegacyFunctionSpecification()  
        function.addParameter('minppbin', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_minppbin():
        """ target no. of particles per time bin"""
        function = LegacyFunctionSpecification()  
        function.addParameter('minppbin', dtype='i', direction=function.OUT)
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
        """ target number of SPH nieghbours"""
        function = LegacyFunctionSpecification()  
        function.addParameter('nsmooth', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function;
    @legacy_function   
    def get_nsmooth():
        """ target number of SPH nieghbours"""
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
        

class Fi(GravitationalDynamics):
    
    def __init__(self, convert_nbody = None):
        
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
       
        
        legacy_interface = BHTreeInterface()
        
        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
        )     
            

class glfi(fi):
    def __init__(self):
        LegacyInterface.__init__(self,name_of_the_worker = 'glworker')
    
    @legacy_function
    def viewer():
        function = LegacyFunctionSpecification()  
        return function
        
