import os
import numpy
from amuse.legacy import *
from amuse.legacy.interface.se import StellarEvolution

from amuse.support.interface import CodeInterface

class MESAInterface(LegacyInterface, LiteratureRefs, StellarEvolution): 
    """
    The software project MESA (Modules for Experiments in Stellar Astrophysics, 
    http://mesa.sourceforge.net/), aims to provide state-of-the-art, robust, 
    and efficient open source modules, usable singly or in combination for a 
    wide range of applications in stellar astrophysics. The AMUSE interface to 
    MESA can create and evolve stars using the MESA/STAR module. If you order a 
    metallicity you haven't used before, starting models will be computed 
    automatically and saved in the `mesa/src/data/star_data/starting_models` 
    directory (please be patient...). All metallicities are supported, even the 
    interesting case of Z=0. The supported stellar mass range is from 
    about 0.1 to 100 Msun.
    
        .. [#] Please acknowledge the use of MESA in your papers. More details
        .. [#] ... on MESA can be found at: http://mesa.sourceforge.net/.
    """
    def __init__(self, **options):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code", **options)
        LiteratureRefs.__init__(self)

    def get_data_directory(self):
        """
        Returns the root name of the directory for the MESA
        application data files.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'mesa', 'input')
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'mesa', 'output')
    
    @property
    def default_path_to_inlist(self):
        return os.path.join(self.get_data_directory(), 'AMUSE_inlist')

    @property
    def default_path_to_MESA_data(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'src', 'data')

    @legacy_function
    def set_MESA_paths():
        """
        Set the paths to the MESA inlist and data directories.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('inlist_path', dtype='string', direction=function.IN,
            description = "Path to the inlist file.")
        function.addParameter('MESA_data_path', dtype='string', direction=function.IN,
            description = "Path to the data directory.")
        function.addParameter('local_data_path', dtype='string', direction=function.IN,
            description = "Path to the data directory.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function
    
    @legacy_function
    def get_maximum_number_of_stars():
        """
        Retrieve the maximum number of stars that can be
        handled by this instance.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('maximum_number_of_stars', dtype='int32', direction=function.OUT,
            description = "The current value of the maximum number of stars")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of was retrieved
        """
        return function
    
    @legacy_function     
    def evolve_to():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('end_time', dtype='float64', direction=function.IN)
        return function
        
    @legacy_function     
    def new_zams_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('status', dtype='int32', direction=function.OUT)
        return function
        
    @legacy_function      
    def get_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('time_step', dtype='float64', direction=function.OUT
            , description="The next timestep for the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
    @legacy_function   
    def get_number_of_zones():
        """
        Retrieve the current number of zones/mesh-cells of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('n_zones', dtype='int32', direction=function.OUT
            , description="The current number of zones/mesh-cells of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_mass_fraction_at_zone():
        """
        Retrieve the mass fraction at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('dq_i', dtype='float64', direction=function.OUT
            , description="The mass fraction at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    @legacy_function   
    def get_temperature_at_zone():
        """
        Retrieve the temperature at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('T_i', dtype='float64', direction=function.OUT
            , description="The temperature at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    @legacy_function   
    def get_density_at_zone():
        """
        Retrieve the density at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('rho_i', dtype='float64', direction=function.OUT
            , description="The density at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    @legacy_function   
    def get_radius_at_zone():
        """
        Retrieve the radius at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('R_i', dtype='float64', direction=function.OUT
            , description="The radius at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    @legacy_function   
    def get_luminosity_at_zone():
        """
        Retrieve the luminosity at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('lum_i', dtype='float64', direction=function.OUT
            , description="The luminosity at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    @legacy_function   
    def get_mu_at_zone():
        """
        Retrieve the mean molecular weight per particle (ions + free electrons)
        at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('mu_i', dtype='float64', direction=function.OUT
            , description="The mean molecular weight at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_number_of_species():
        """
        Retrieve the current number of chemical abundance variables per zone of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('n_species', dtype='int32', direction=function.OUT
            , description="The current number of chemical abundance variables per zone of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_name_of_species():
        """
        Retrieve the name of the chemical abundance variable of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to get the name of")
        function.addParameter('species_name', dtype='string', direction=function.OUT
            , description="The name of the chemical abundance variable of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_id_of_species():
        """
        Retrieve the chem_ID of the chemical abundance variable of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to get the name of")
        function.addParameter('species_id', dtype='int32', direction=function.OUT
            , description="The chem_ID of the chemical abundance variable of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_mass_of_species():
        """
        Retrieve the mass number of the chemical abundance variable of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to get the mass number of")
        function.addParameter('species_mass', dtype='float64', direction=function.OUT
            , description="The mass number of the chemical abundance variable of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_mass_fraction_of_species_at_zone():
        """
        Retrieve the fractional chemical abundance variable at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('Xj_i', dtype='float64', direction=function.OUT
            , description="The fractional chemical abundance variable at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function
    def get_max_age_stop_condition():
        """
        Retrieve the current maximum age stop condition of this instance (in years).
        Evolution will stop once the star has reached this maximum age.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('max_age_stop_condition', dtype='float64', direction=function.OUT
            , description="The current maximum age stop condition of this instance (in years).")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_max_age_stop_condition():
        """
        Set the new maximum age stop condition of this instance (in years).
        Evolution will stop once the star has reached this maximum age.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('max_age_stop_condition', dtype='float64', direction=function.IN
            , description="The new maximum age stop condition of this instance (in years).")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_min_timestep_stop_condition():
        """
        Retrieve the current minimum timestep stop condition of this instance (in years).
        Evolution will stop if the timestep required by the solver in order to converge
        has decreased below this minimum timestep.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('min_timestep_stop_condition', dtype='float64', direction=function.OUT
            , description="The current minimum timestep stop condition of this instance (in years).")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_min_timestep_stop_condition():
        """
        Set the new minimum timestep stop condition of this instance (in years).
        Evolution will stop if the timestep required by the solver in order to converge
        has decreased below this minimum timestep.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('min_timestep_stop_condition', dtype='float64', direction=function.IN
            , description="The new minimum timestep stop condition of this instance (in years).")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_max_iter_stop_condition():
        """
        Retrieve the current maximum number of iterations of this instance. (Negative means no maximum)
        Evolution will stop after this number of iterations.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('max_iter_stop_condition', dtype='int32', direction=function.OUT
            , description="The current maximum number of iterations of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_max_iter_stop_condition():
        """
        Set the new maximum number of iterations of this instance. (Negative means no maximum)
        Evolution will stop after this number of iterations.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('max_iter_stop_condition', dtype='int32', direction=function.IN
            , description="The new maximum number of iterations of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_mixing_length_ratio():
        """
        Retrieve the current value of the mixing length ratio.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('mixing_length_ratio', dtype='float64', direction=function.OUT
            , description="The current value of the mixing length ratio.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_mixing_length_ratio():
        """
        Set the value of the mixing length ratio.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('mixing_length_ratio', dtype='float64', direction=function.IN
            , description="The new value of the mixing length ratio.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_semi_convection_efficiency():
        """
        Retrieve the current value of the efficiency of semi-convection,
        after Heger, Langer, & Woosley 2000 (ApJ), which goes back to 
        Langer, Sugimoto & Fricke 1983 (A&A).
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('semi_convection_efficiency', dtype='float64', direction=function.OUT
            , description="The current value of the efficiency of semi-convection.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_semi_convection_efficiency():
        """
        Set the value of the efficiency of semi-convection,
        after Heger, Langer, & Woosley 2000 (ApJ), which goes back to 
        Langer, Sugimoto & Fricke 1983 (A&A).
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('semi_convection_efficiency', dtype='float64', direction=function.IN
            , description="The new value of the efficiency of semi-convection.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
    
    @legacy_function
    def get_RGB_wind_scheme():
        """
        Retrieve the current wind (mass loss) scheme for RGB stars:
        No automatic wind (0)
        Reimers (1): e.g. see: Baschek, Kegel, Traving (eds), Springer, Berlin, 1975, p. 229.
        Blocker (2): T. Blocker, A&A 297, 727-738 (1995)
        de Jager (3): de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259
        Dutch (4): Glebbeek et al 2009, Vink et al 2001, Nugis & Lamers 2000, de Jager 1990
        Mattsson (5)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('RGB_wind_scheme', dtype='int32', direction=function.OUT
            , description="The current wind (mass loss) scheme for RGB stars of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_RGB_wind_scheme():
        """
        Set the new wind (mass loss) scheme for RGB stars:
        No automatic wind (0)
        Reimers (1): e.g. see: Baschek, Kegel, Traving (eds), Springer, Berlin, 1975, p. 229.
        Blocker (2): T. Blocker, A&A 297, 727-738 (1995)
        de Jager (3): de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259
        Dutch (4): Glebbeek et al 2009, Vink et al 2001, Nugis & Lamers 2000, de Jager 1990
        Mattsson (5)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('RGB_wind_scheme', dtype='int32', direction=function.IN
            , description="The new wind (mass loss) scheme for RGB stars of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
    
    @legacy_function
    def get_RGB_wind_efficiency():
        """
        Retrieve the current mass loss efficiency for RGB stars.
        Exact implementation depends on RGB_wind_scheme.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('RGB_wind_efficiency', dtype='float64', direction=function.OUT
            , description="The current value of the mass loss efficiency for RGB stars.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_RGB_wind_efficiency():
        """
        Set the value of the mass loss efficiency for RGB stars.
        Exact implementation depends on RGB_wind_scheme.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('RGB_wind_efficiency', dtype='float64', direction=function.IN
            , description="The new value of the mass loss efficiency for RGB stars.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
    
    @legacy_function
    def get_AGB_wind_scheme():
        """
        Retrieve the current wind (mass loss) scheme for AGB stars:
        No automatic wind (0)
        Reimers (1): e.g. see: Baschek, Kegel, Traving (eds), Springer, Berlin, 1975, p. 229.
        Blocker (2): T. Blocker, A&A 297, 727-738 (1995)
        de Jager (3): de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259
        Dutch (4): Glebbeek et al 2009, Vink et al 2001, Nugis & Lamers 2000, de Jager 1990
        Mattsson (5)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('AGB_wind_scheme', dtype='int32', direction=function.OUT
            , description="The current wind (mass loss) scheme for AGB stars of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_AGB_wind_scheme():
        """
        Set the new wind (mass loss) scheme for AGB stars:
        No automatic wind (0)
        Reimers (1): e.g. see: Baschek, Kegel, Traving (eds), Springer, Berlin, 1975, p. 229.
        Blocker (2): T. Blocker, A&A 297, 727-738 (1995)
        de Jager (3): de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259
        Dutch (4): Glebbeek et al 2009, Vink et al 2001, Nugis & Lamers 2000, de Jager 1990
        Mattsson (5)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('AGB_wind_scheme', dtype='int32', direction=function.IN
            , description="The new wind (mass loss) scheme for AGB stars of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
    
    @legacy_function
    def get_AGB_wind_efficiency():
        """
        Retrieve the current mass loss efficiency for AGB stars.
        Exact implementation depends on AGB_wind_scheme.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('AGB_wind_efficiency', dtype='float64', direction=function.OUT
            , description="The current value of the mass loss efficiency for AGB stars.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_AGB_wind_efficiency():
        """
        Set the value of the mass loss efficiency for AGB stars.
        Exact implementation depends on AGB_wind_scheme.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('AGB_wind_efficiency', dtype='float64', direction=function.IN
            , description="The new value of the mass loss efficiency for AGB stars.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

class MESA(CodeInterface):
    
    def __init__(self, **options):
        CodeInterface.__init__(self, MESAInterface(), **options)
        self.set_MESA_paths(self.default_path_to_inlist, 
            self.default_path_to_MESA_data, self.get_data_directory())
        self.parameters.set_defaults()
        
    
    def define_parameters(self, object):
        
        object.add_method_parameter(
            "get_metallicity",
            "set_metallicity",
            "metallicity", 
            "Metallicity of all stars", 
            units.none, 
            0.02 | units.none
        )
        
        object.add_method_parameter(
            "get_max_age_stop_condition",
            "set_max_age_stop_condition",
            "max_age_stop_condition", 
            "The maximum age stop condition of this instance.",
            units.yr, 
            1.0e12 | units.yr
        )
        
        object.add_method_parameter(
            "get_min_timestep_stop_condition",
            "set_min_timestep_stop_condition",
            "min_timestep_stop_condition", 
            "The minimum timestep stop condition of this instance.",
            units.yr, 
            1.0e6 | units.s
        )
        
        object.add_method_parameter(
            "get_max_iter_stop_condition",
            "set_max_iter_stop_condition",
            "max_iter_stop_condition", 
            "The maximum number of iterations of this instance. (Negative means no maximum)",
            units.none, 
            -1111 | units.none
        )
        
        object.add_method_parameter(
            "get_mixing_length_ratio",
            "set_mixing_length_ratio",
            "mixing_length_ratio", 
            "The mixing-length ratio (alpha).",
            units.none, 
            2.0 | units.none
        )
        
        object.add_method_parameter(
            "get_semi_convection_efficiency",
            "set_semi_convection_efficiency",
            "semi_convection_efficiency", 
            "The efficiency of semi-convection, after Heger, Langer, & Woosley 2000 (ApJ), "
               "which goes back to Langer, Sugimoto & Fricke 1983 (A&A).",
            units.none, 
            0.0 | units.none
        )
        
        object.add_method_parameter(
            "get_RGB_wind_scheme",
            "set_RGB_wind_scheme",
            "RGB_wind_scheme", 
            "The mass loss scheme for RGB stars: none (0), Reimers (1), "
                "Blocker (2), de Jager (3), Dutch (4), Mattsson (5)",
            units.none, 
            0 | units.none
        )
        
        object.add_method_parameter(
            "get_AGB_wind_scheme",
            "set_AGB_wind_scheme",
            "AGB_wind_scheme", 
            "The mass loss scheme for AGB stars: none (0), Reimers (1), "
                "Blocker (2), de Jager (3), Dutch (4), Mattsson (5)",
            units.none, 
            0 | units.none
        )
        
        object.add_method_parameter(
            "get_RGB_wind_efficiency",
            "set_RGB_wind_efficiency",
            "RGB_wind_efficiency", 
            "The mass loss efficiency for RGB stars. Exact implementation depends on RGB_wind_scheme.",
            units.none, 
            0.0 | units.none
        )
        
        object.add_method_parameter(
            "get_AGB_wind_efficiency",
            "set_AGB_wind_efficiency",
            "AGB_wind_efficiency", 
            "The mass loss efficiency for AGB stars. Exact implementation depends on AGB_wind_scheme.",
            units.none, 
            0.0 | units.none
        )
        
        
        
    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_star')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_star')
        
        object.add_getter('particles', 'get_radius', names = ('radius',))
        object.add_getter('particles', 'get_stellar_type', names = ('stellar_type',))
        object.add_getter('particles', 'get_mass', names = ('mass',))
        object.add_getter('particles', 'get_age', names = ('age',))
        object.add_getter('particles', 'get_time_step', names = ('time_step',))
        object.add_getter('particles', 'get_luminosity',names = ('luminosity',))
        object.add_getter('particles', 'get_temperature',names = ('temperature',))
        
        object.add_method('particles', 'get_number_of_zones')
        object.add_method('particles', 'get_mass_profile')
        object.add_method('particles', 'get_density_profile')
        object.add_method('particles', 'get_radius_profile')
        object.add_method('particles', 'get_temperature_profile')
        object.add_method('particles', 'get_luminosity_profile')
        object.add_method('particles', 'get_mu_profile')
        object.add_method('particles', 'get_number_of_species')
        object.add_method('particles', 'get_names_of_species')
        object.add_method('particles', 'get_IDs_of_species')
        object.add_method('particles', 'get_masses_of_species')
        object.add_method('particles', 'get_chemical_abundance_profiles')
        object.add_method('particles', 'evolve', 'evolve_one_step')
    
    def define_errorcodes(self, object):
        object.add_errorcode(-1, 'Something went wrong...')
        object.add_errorcode(-11, 'Evolve terminated: Unspecified stop condition reached.')
        object.add_errorcode(-12, 'Evolve terminated: Maximum age reached.')
        object.add_errorcode(-13, 'Evolve terminated: Maximum number of iterations reached.')
        object.add_errorcode(-14, 'Evolve terminated: Maximum number of backups reached.')
    
    def define_methods(self, object):
        
        object.add_method(
            "evolve",
            (object.INDEX,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "new_particle",
            (units.MSun),
            (object.INDEX, object.ERROR_CODE)
        )
        object.add_method(
            "delete_star",
            (object.INDEX,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_mass",
            (object.INDEX,),
            (units.MSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_radius",
            (object.INDEX,),
            (units.RSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_stellar_type",
            (object.INDEX,),
            (units.stellar_type, object.ERROR_CODE,)
        )
        object.add_method(
            "get_age", 
            (object.INDEX,), 
            (units.yr, object.ERROR_CODE,)
        )
        object.add_method(
            "get_luminosity", 
            (object.INDEX,), 
            (units.LSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_temperature", 
            (object.INDEX,), 
            (units.K, object.ERROR_CODE,)
        )
        object.add_method(
            "get_time_step", 
            (object.INDEX,), 
            (units.yr, object.ERROR_CODE,)
        )
        object.add_method(
            "get_number_of_zones", 
            (object.INDEX,), 
            (units.none, object.ERROR_CODE,)
        )
        object.add_method(
            "get_mass_fraction_at_zone", 
            (object.INDEX,units.none,), 
            (units.none, object.ERROR_CODE,)
        )
        object.add_method(
            "get_temperature_at_zone", 
            (object.INDEX,units.none,), 
            (units.K, object.ERROR_CODE,)
        )
        object.add_method(
            "get_density_at_zone", 
            (object.INDEX,units.none,), 
            (units.g/units.cm**3, object.ERROR_CODE,)
        )
        object.add_method(
            "get_radius_at_zone", 
            (object.INDEX,units.none,), 
            (units.cm, object.ERROR_CODE,)
        )
        object.add_method(
            "get_luminosity_at_zone", 
            (object.INDEX,units.none,), 
            (units.erg/units.s, object.ERROR_CODE,)
        )
        object.add_method(
            "get_mu_at_zone", 
            (object.INDEX,units.none,), 
            (units.amu, object.ERROR_CODE,)
        )
        object.add_method(
            "get_number_of_species", 
            (object.INDEX,), 
            (units.none, object.ERROR_CODE,)
        )
        object.add_method(
            "get_name_of_species", 
            (object.INDEX,units.none,), 
            (units.string, object.ERROR_CODE,)
        )
        object.add_method(
            "get_id_of_species", 
            (object.INDEX,units.none,), 
            (units.none, object.ERROR_CODE,)
        )
        object.add_method(
            "get_mass_of_species", 
            (object.INDEX,units.none,), 
            (units.amu, object.ERROR_CODE,)
        )
        object.add_method(
            "get_mass_fraction_of_species_at_zone", 
            (object.INDEX,units.none,units.none,), 
            (units.none, object.ERROR_CODE,)
        )
        
    
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
        self.initialize_code()
        
    def initialize_module_with_current_parameters(self):
        self.initialize_code()
        
    
        
    def evolve_model(self, end_time = None):
        if end_time is None:
            result = self.particles.evolve_one_step()
            return result
                   
        for particle in self.particles:
            while particle.age < end_time:
                particle.evolve_one_step()
    
    def get_mass_profile(self, indices_of_the_stars, number_of_zones = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise exceptions.LegacyException("Querying mass profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_mass_fraction_at_zone([indices_of_the_stars]*number_of_zones, range(1,number_of_zones+1) | units.none)
    
    def get_density_profile(self, indices_of_the_stars, number_of_zones = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying density profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_density_at_zone([indices_of_the_stars]*number_of_zones, range(1,number_of_zones+1) | units.none)
    
    def get_radius_profile(self, indices_of_the_stars, number_of_zones = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying radius profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_radius_at_zone([indices_of_the_stars]*number_of_zones, range(1,number_of_zones+1) | units.none)
    
    def get_temperature_profile(self, indices_of_the_stars, number_of_zones = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying temperature profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_temperature_at_zone([indices_of_the_stars]*number_of_zones, range(1,number_of_zones+1) | units.none)
    
    def get_luminosity_profile(self, indices_of_the_stars, number_of_zones = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying luminosity profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_luminosity_at_zone([indices_of_the_stars]*number_of_zones, range(1,number_of_zones+1) | units.none)
    
    def get_mu_profile(self, indices_of_the_stars, number_of_zones = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying mean-molecular-weight profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_mu_at_zone([indices_of_the_stars]*number_of_zones, range(1,number_of_zones+1) | units.none)
    
    
    def get_names_of_species(self, indices_of_the_stars, number_of_species = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying chemical abundance names of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        return list(self.get_name_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1) | units.none
        ).value_in(units.string))
    
    def get_IDs_of_species(self, indices_of_the_stars, number_of_species = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying chemical abundance IDs of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        return list(self.get_id_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1) | units.none
        ).value_in(units.none))
    
    def get_masses_of_species(self, indices_of_the_stars, number_of_species = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying chemical abundance mass numbers of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        return self.get_mass_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1) | units.none
        )
    
    
    def get_chemical_abundance_profiles(self, indices_of_the_stars, number_of_zones = None, number_of_species = None):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise Exception("Querying chemical abundance profiles of more than one particle at a time is not supported.")
            indices_of_the_stars = indices_of_the_stars[0]
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        grid = numpy.indices((number_of_species, number_of_zones)) + 1
        return self.get_mass_fraction_of_species_at_zone(
            [indices_of_the_stars] * number_of_zones * number_of_species, 
            units.none.new_quantity(grid[0].flatten()), 
            units.none.new_quantity(grid[1].flatten())
        ).reshape((number_of_species, number_of_zones))
    
