import os
import numpy
from operator import itemgetter

from amuse.community import *
from amuse.community.interface.se import StellarEvolution, StellarEvolutionInterface, \
    InternalStellarStructure, InternalStellarStructureInterface

from amuse.units.quantities import VectorQuantity
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.options import option

class MESAInterface(CodeInterface, LiteratureReferencesMixIn, StellarEvolutionInterface, 
        InternalStellarStructureInterface, CodeWithDataDirectories): 
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
    
    References:
        .. [#] Paxton, Bildsten, Dotter, Herwig, Lesaffre & Timmes 2010, ApJS submitted, arXiv:1009.1622
        .. [#] http://mesa.sourceforge.net/
    """
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="mesa_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    @property
    def default_path_to_inlist(self):
        return os.path.join(self.get_data_directory(), 'AMUSE_inlist')

    @option(type="string", sections=('data'))
    def default_path_to_MESA_data(self):
        return os.path.join(self.amuse_root_directory, 'src', 'amuse', 'community', 'mesa', 'src', 'mesa', 'data')
    
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
    def new_zams_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('status', dtype='int32', direction=function.OUT)
        return function
        
    @legacy_function   
    def new_pre_ms_particle():
        """
        Define a new pre-main-sequence star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('mass', dtype='float64', direction=function.IN
            , description="The initial mass of the star")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('time_step', dtype='float64', direction=function.IN
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
    def get_core_mass():
        """
        Retrieve the current core mass of the star, where hydrogen abundance is <= h1_boundary_limit
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('core_mass', dtype='float64', direction=function.OUT
            , description="The current core mass of the star, where hydrogen abundance is <= h1_boundary_limit")
        function.result_type = 'int32'
        return function
    
    @legacy_function   
    def get_mass_loss_rate():
        """
        Retrieve the current mass loss rate of the star. (positive for winds, negative for accretion)
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('mass_loss_rate', dtype='float64', direction=function.OUT
            , description="The current mass loss rate of the star. (positive for winds, negative for accretion)")
        function.result_type = 'int32'
        return function
    
    @legacy_function   
    def get_manual_mass_transfer_rate():
        """
        Retrieve the current user-specified mass transfer rate of the star. (negative for winds, positive for accretion)
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('mass_change', dtype='float64', direction=function.OUT
            , description="The current user-specified mass transfer rate of the star. (negative for winds, positive for accretion)")
        function.result_type = 'int32'
        return function
    
    @legacy_function   
    def set_manual_mass_transfer_rate():
        """
        Set a new user-specified mass transfer rate of the star. (negative for winds, positive for accretion)
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('mass_change', dtype='float64', direction=function.IN
            , description="The new user-specified mass transfer rate of the star. (negative for winds, positive for accretion)")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_accrete_same_as_surface():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('accrete_same_as_surface_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_accrete_same_as_surface():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('accrete_same_as_surface_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_accrete_composition_non_metals():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('h1', dtype='float64', direction=function.OUT)
        function.addParameter('h2', dtype='float64', direction=function.OUT)
        function.addParameter('he3', dtype='float64', direction=function.OUT)
        function.addParameter('he4', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_accrete_composition_non_metals():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('h1', dtype='float64', direction=function.IN)
        function.addParameter('h2', dtype='float64', direction=function.IN)
        function.addParameter('he3', dtype='float64', direction=function.IN)
        function.addParameter('he4', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_accrete_composition_metals_identifier():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('accrete_composition_metals_identifier', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_accrete_composition_metals_identifier():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('accrete_composition_metals_identifier', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_accrete_composition_metals():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('li', dtype='float64', direction=function.OUT)
        function.addParameter('be', dtype='float64', direction=function.OUT)
        function.addParameter('b', dtype='float64', direction=function.OUT)
        function.addParameter('c', dtype='float64', direction=function.OUT)
        function.addParameter('n', dtype='float64', direction=function.OUT)
        function.addParameter('o', dtype='float64', direction=function.OUT)
        function.addParameter('f', dtype='float64', direction=function.OUT)
        function.addParameter('ne', dtype='float64', direction=function.OUT)
        function.addParameter('na', dtype='float64', direction=function.OUT)
        function.addParameter('mg', dtype='float64', direction=function.OUT)
        function.addParameter('al', dtype='float64', direction=function.OUT)
        function.addParameter('si', dtype='float64', direction=function.OUT)
        function.addParameter('p', dtype='float64', direction=function.OUT)
        function.addParameter('s', dtype='float64', direction=function.OUT)
        function.addParameter('cl', dtype='float64', direction=function.OUT)
        function.addParameter('ar', dtype='float64', direction=function.OUT)
        function.addParameter('k', dtype='float64', direction=function.OUT)
        function.addParameter('ca', dtype='float64', direction=function.OUT)
        function.addParameter('sc', dtype='float64', direction=function.OUT)
        function.addParameter('ti', dtype='float64', direction=function.OUT)
        function.addParameter('v', dtype='float64', direction=function.OUT)
        function.addParameter('cr', dtype='float64', direction=function.OUT)
        function.addParameter('mn', dtype='float64', direction=function.OUT)
        function.addParameter('fe', dtype='float64', direction=function.OUT)
        function.addParameter('co', dtype='float64', direction=function.OUT)
        function.addParameter('ni', dtype='float64', direction=function.OUT)
        function.addParameter('cu', dtype='float64', direction=function.OUT)
        function.addParameter('zn', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_accrete_composition_metals():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN)
        function.addParameter('li', dtype='float64', direction=function.IN)
        function.addParameter('be', dtype='float64', direction=function.IN)
        function.addParameter('b', dtype='float64', direction=function.IN)
        function.addParameter('c', dtype='float64', direction=function.IN)
        function.addParameter('n', dtype='float64', direction=function.IN)
        function.addParameter('o', dtype='float64', direction=function.IN)
        function.addParameter('f', dtype='float64', direction=function.IN)
        function.addParameter('ne', dtype='float64', direction=function.IN)
        function.addParameter('na', dtype='float64', direction=function.IN)
        function.addParameter('mg', dtype='float64', direction=function.IN)
        function.addParameter('al', dtype='float64', direction=function.IN)
        function.addParameter('si', dtype='float64', direction=function.IN)
        function.addParameter('p', dtype='float64', direction=function.IN)
        function.addParameter('s', dtype='float64', direction=function.IN)
        function.addParameter('cl', dtype='float64', direction=function.IN)
        function.addParameter('ar', dtype='float64', direction=function.IN)
        function.addParameter('k', dtype='float64', direction=function.IN)
        function.addParameter('ca', dtype='float64', direction=function.IN)
        function.addParameter('sc', dtype='float64', direction=function.IN)
        function.addParameter('ti', dtype='float64', direction=function.IN)
        function.addParameter('v', dtype='float64', direction=function.IN)
        function.addParameter('cr', dtype='float64', direction=function.IN)
        function.addParameter('mn', dtype='float64', direction=function.IN)
        function.addParameter('fe', dtype='float64', direction=function.IN)
        function.addParameter('co', dtype='float64', direction=function.IN)
        function.addParameter('ni', dtype='float64', direction=function.IN)
        function.addParameter('cu', dtype='float64', direction=function.IN)
        function.addParameter('zn', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_number_of_backups_in_a_row():
        """
        Retrieve the number_of_backups_in_a_row of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of number_of_backups_in_a_row")
        function.addParameter('n_backup', dtype='int32', direction=function.OUT
            , description="The current number_of_backups_in_a_row of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The number_of_backups_in_a_row was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    
    @legacy_function
    def reset_number_of_backups_in_a_row():
        """
        Reset number_of_backups_in_a_row of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to reset the value of number_of_backups_in_a_row")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The number_of_backups_in_a_row was reset.
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
    def set_mass_fraction_at_zone():
        """
        Set the mass fraction at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to set the value of")
        function.addParameter('dq_i', dtype='float64', direction=function.IN
            , description="The mass fraction at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
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
    def set_luminosity_at_zone():
        """
        Set the luminosity at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to set the value of")
        function.addParameter('lum_i', dtype='float64', direction=function.IN
            , description="The luminosity at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A zone with the given index was not found.
        """
        return function
    
    @legacy_function
    def get_entropy_at_zone():
        """
        Retrieve the entropy at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('S_i', dtype='float64', direction=function.OUT
            , description="The specific entropy at the specified zone/mesh-cell of the star.")
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
    def get_thermal_energy_at_zone():
        """
        Retrieve the entropy at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('E_i', dtype='float64', direction=function.OUT
            , description="The specific thermal energy at the specified zone/mesh-cell of the star.")
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
    def get_brunt_vaisala_frequency_squared_at_zone():
        """
        Retrieve the Brunt-Vaisala frequency squared at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN, unit=INDEX)
        function.addParameter('zone', dtype='int32', direction=function.IN, unit=NO_UNIT)
        function.addParameter('brunt_N2', dtype='float64', direction=function.OUT, unit=units.s**-2)
        function.result_type = 'int32'
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
    def erase_memory():
        """
        Erase memory of the star, i.e. copy the current structure over the memory of 
        the structure of the previous steps. Useful after setting the stucture of 
        the star, to prevent backup steps to undo changes
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.result_type = 'int32'
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
    def get_convective_overshoot_parameter():
        function = LegacyFunctionSpecification()  
        function.addParameter('convective_overshoot_parameter', dtype='float64', direction=function.OUT,
            description="The current value of the convective overshoot parameter.")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_convective_overshoot_parameter():
        function = LegacyFunctionSpecification()  
        function.addParameter('convective_overshoot_parameter', dtype='float64', direction=function.IN,
            description="The new value of the convective overshoot parameter.")
        function.result_type = 'int32'
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
    def get_reimers_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('reimers_wind_efficiency', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_reimers_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('reimers_wind_efficiency', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_blocker_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('blocker_wind_efficiency', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_blocker_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('blocker_wind_efficiency', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_de_jager_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('de_jager_wind_efficiency', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_de_jager_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('de_jager_wind_efficiency', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_dutch_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('dutch_wind_efficiency', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_dutch_wind_efficiency():
        function = LegacyFunctionSpecification()  
        function.addParameter('dutch_wind_efficiency', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_stabilize_new_stellar_model_flag():
        function = LegacyFunctionSpecification()  
        function.addParameter('stabilize_new_stellar_model_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_stabilize_new_stellar_model_flag():
        function = LegacyFunctionSpecification()  
        function.addParameter('stabilize_new_stellar_model_flag', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function   
    def new_stellar_model():
        """
        Define a new star model in the code. The star needs to be finalized 
        before it can evolve, see 'finalize_stellar_model'.
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for par in ['d_mass', 'radius', 'rho', 'temperature', 'luminosity', 
                'X_H', 'X_He', 'X_C', 'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']:
            function.addParameter(par, dtype='float64', direction=function.IN)
        function.addParameter('n', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function   
    def finalize_stellar_model():
        """
        Finalize the new star model defined by 'new_stellar_model'.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.OUT, description = "The new index for the star. "
            "This index can be used to refer to this star in other functions")
        function.addParameter('age_tag', dtype='float64', direction=function.IN, 
            description = "The initial age of the star")
        function.result_type = 'int32'
        return function

class MESA(StellarEvolution, InternalStellarStructure):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, MESAInterface(**options), **options)
        
        output_dir = self.get_output_directory()
        if not self.channel_type == 'distributed':
            self.ensure_data_directory_exists(os.path.join(output_dir, 'star_data', 'starting_models'))
        
        self.set_MESA_paths(
            self.default_path_to_inlist, 
            self.default_path_to_MESA_data, 
            output_dir
        )
        self.model_time = 0.0 | units.yr
        
    
    def define_parameters(self, handler):
        
        handler.add_method_parameter(
            "get_metallicity",
            "set_metallicity",
            "metallicity", 
            "Metallicity of all stars", 
            default_value = 0.02
        )
        
        handler.add_method_parameter(
            "get_max_age_stop_condition",
            "set_max_age_stop_condition",
            "max_age_stop_condition", 
            "The maximum age stop condition of this instance.",
            default_value = 1.0e36 | units.yr
        )
        
        handler.add_method_parameter(
            "get_min_timestep_stop_condition",
            "set_min_timestep_stop_condition",
            "min_timestep_stop_condition", 
            "The minimum timestep stop condition of this instance.",
            default_value = 1.0e-6 | units.s
        )
        
        handler.add_method_parameter(
            "get_max_iter_stop_condition",
            "set_max_iter_stop_condition",
            "max_iter_stop_condition", 
            "The maximum number of iterations of this instance. (Negative means no maximum)",
            default_value = -1111
        )
        
        handler.add_method_parameter(
            "get_convective_overshoot_parameter",
            "set_convective_overshoot_parameter",
            "herwig_convective_overshoot_parameter", 
            "The convective overshoot parameter (Herwig 2000), f=0.016 is argued to be a reasonable value.",
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_mixing_length_ratio",
            "set_mixing_length_ratio",
            "mixing_length_ratio", 
            "The mixing-length ratio (alpha).",
            default_value = 2.0
        )
        
        handler.add_method_parameter(
            "get_semi_convection_efficiency",
            "set_semi_convection_efficiency",
            "semi_convection_efficiency", 
            "The efficiency of semi-convection, after Heger, Langer, & Woosley 2000 (ApJ), "
               "which goes back to Langer, Sugimoto & Fricke 1983 (A&A).",
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_RGB_wind_scheme",
            "set_RGB_wind_scheme",
            "RGB_wind_scheme", 
            "The mass loss scheme for RGB stars: none (0), Reimers (1), "
                "Blocker (2), de Jager (3), Dutch (4), Mattsson (5)",
            default_value = 1
        )
        
        handler.add_method_parameter(
            "get_AGB_wind_scheme",
            "set_AGB_wind_scheme",
            "AGB_wind_scheme", 
            "The mass loss scheme for AGB stars: none (0), Reimers (1), "
                "Blocker (2), de Jager (3), Dutch (4), Mattsson (5)",
            default_value = 1
        )
        
        handler.add_method_parameter(
            "get_reimers_wind_efficiency",
            "set_reimers_wind_efficiency",
            "reimers_wind_efficiency", 
            "The Reimers mass loss efficiency. Only used if (RGB/AGB_wind_scheme == 1).",
            default_value = 0.5
        )
        handler.add_method_parameter(
            "get_blocker_wind_efficiency",
            "set_blocker_wind_efficiency",
            "blocker_wind_efficiency", 
            "The Blocker mass loss efficiency. Only used if (RGB/AGB_wind_scheme == 2).",
            default_value = 0.1
        )
        handler.add_method_parameter(
            "get_de_jager_wind_efficiency",
            "set_de_jager_wind_efficiency",
            "de_jager_wind_efficiency", 
            "The de Jager mass loss efficiency. Only used if (RGB/AGB_wind_scheme == 3).",
            default_value = 0.8
        )
        handler.add_method_parameter(
            "get_dutch_wind_efficiency",
            "set_dutch_wind_efficiency",
            "dutch_wind_efficiency", 
            "The Dutch mass loss efficiency. Only used if (RGB/AGB_wind_scheme == 4).",
            default_value = 0.8
        )
        handler.add_boolean_parameter(
            "get_stabilize_new_stellar_model_flag",
            "set_stabilize_new_stellar_model_flag",
            "stabilize_new_stellar_model_flag",
            "Flag specifying whether to stabilize any loaded stellar models first.",
            default_value = True
        )
        
        
        
    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['native_stars', 'imported_stars', 'pre_ms_stars'], 
            index_to_default_set = 0)
        
        handler.define_set('imported_stars', 'index_of_the_star')
        handler.set_new('imported_stars', 'finalize_stellar_model')
        handler.set_delete('imported_stars', 'delete_star')
        
        handler.define_set('native_stars', 'index_of_the_star')
        handler.set_new('native_stars', 'new_particle')
        handler.set_delete('native_stars', 'delete_star')
        
        handler.define_set('pre_ms_stars', 'index_of_the_star')
        handler.set_new('pre_ms_stars', 'new_pre_ms_particle')
        handler.set_delete('pre_ms_stars', 'delete_star')
        
        for particle_set_name in ['native_stars', 'imported_stars', 'pre_ms_stars']:
            handler.add_getter(particle_set_name, 'get_radius', names = ('radius',))
            handler.add_getter(particle_set_name, 'get_stellar_type', names = ('stellar_type',))
            handler.add_getter(particle_set_name, 'get_mass', names = ('mass',))
            handler.add_setter(particle_set_name, 'set_mass', names = ('mass',))
            handler.add_getter(particle_set_name, 'get_core_mass', names = ('core_mass',))
            handler.add_getter(particle_set_name, 'get_mass_loss_rate', names = ('wind',))
            handler.add_getter(particle_set_name, 'get_age', names = ('age',))
            handler.add_getter(particle_set_name, 'get_time_step', names = ('time_step',))
            handler.add_setter(particle_set_name, 'set_time_step', names = ('time_step',))
            handler.add_getter(particle_set_name, 'get_luminosity', names = ('luminosity',))
            handler.add_getter(particle_set_name, 'get_temperature', names = ('temperature',))
            
            handler.add_getter(particle_set_name, 'get_manual_mass_transfer_rate', names = ('mass_change',))
            handler.add_setter(particle_set_name, 'set_manual_mass_transfer_rate', names = ('mass_change',))
            
            handler.add_method(particle_set_name, 'get_accrete_same_as_surface')
            handler.add_method(particle_set_name, 'set_accrete_same_as_surface')
            handler.add_method(particle_set_name, 'get_accrete_composition_non_metals')
            handler.add_method(particle_set_name, 'set_accrete_composition_non_metals')
            handler.add_method(particle_set_name, 'get_accrete_composition_metals_identifier')
            handler.add_method(particle_set_name, 'set_accrete_composition_metals_identifier')
            handler.add_method(particle_set_name, 'get_accrete_composition_metals')
            handler.add_method(particle_set_name, 'set_accrete_composition_metals')
            
            handler.add_method(particle_set_name, 'evolve_one_step')
            handler.add_method(particle_set_name, 'evolve_for')
            InternalStellarStructure.define_particle_sets(
                self, 
                handler, 
                set_name = particle_set_name
            )
            handler.add_method(particle_set_name, 'get_mass_profile')
            handler.add_method(particle_set_name, 'set_mass_profile')
            handler.add_method(particle_set_name, 'get_cumulative_mass_profile')
            handler.add_method(particle_set_name, 'get_luminosity_profile')
            handler.add_method(particle_set_name, 'set_luminosity_profile')
            handler.add_method(particle_set_name, 'get_entropy_profile')
            handler.add_method(particle_set_name, 'get_thermal_energy_profile')
            handler.add_method(particle_set_name, 'get_brunt_vaisala_frequency_squared_profile')
            handler.add_method(particle_set_name, 'get_IDs_of_species')
            handler.add_method(particle_set_name, 'get_masses_of_species')
            handler.add_method(particle_set_name, 'get_number_of_backups_in_a_row')
            handler.add_method(particle_set_name, 'reset_number_of_backups_in_a_row')
            
    def define_state(self, handler):
        StellarEvolution.define_state(self, handler)
        handler.add_method('EDIT', 'new_pre_ms_particle')
        handler.add_method('UPDATE', 'new_pre_ms_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_pre_ms_particle', False)
        handler.add_method('EDIT', 'finalize_stellar_model')
        handler.add_method('UPDATE', 'finalize_stellar_model')
        handler.add_transition('RUN', 'UPDATE', 'finalize_stellar_model', False)
    
    def define_errorcodes(self, handler):
        InternalStellarStructure.define_errorcodes(self, handler)
        handler.add_errorcode(-1, 'Something went wrong...')
        handler.add_errorcode(-4, 'Not implemented.')
        handler.add_errorcode(-11, 'Evolve terminated: Unspecified stop condition reached.')
        handler.add_errorcode(-12, 'Evolve terminated: Maximum age reached.')
        handler.add_errorcode(-13, 'Evolve terminated: Maximum number of iterations reached.')
        handler.add_errorcode(-14, 'Evolve terminated: Maximum number of backups reached.')
        handler.add_errorcode(-15, 'Evolve terminated: Minimum timestep limit reached.')
    
    def define_methods(self, handler):
        InternalStellarStructure.define_methods(self, handler)
        StellarEvolution.define_methods(self, handler)
        handler.add_method(
            "new_pre_ms_particle",
            (units.MSun),
            (handler.INDEX, handler.ERROR_CODE)
        )
        handler.add_method(
            "set_time_step", 
            (handler.INDEX, units.yr), 
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_core_mass",
            (handler.INDEX,),
            (units.MSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_mass_loss_rate",
            (handler.INDEX,),
            (units.g / units.s, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_manual_mass_transfer_rate",
            (handler.INDEX,),
            (units.MSun / units.yr, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_manual_mass_transfer_rate",
            (handler.INDEX, units.MSun / units.yr),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_number_of_backups_in_a_row", 
            (handler.INDEX,), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "reset_number_of_backups_in_a_row", 
            (handler.INDEX,), 
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_mass_fraction_at_zone", 
            (handler.INDEX,handler.NO_UNIT,), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_mass_fraction_at_zone", 
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT,), 
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_luminosity_at_zone", 
            (handler.INDEX,handler.NO_UNIT,), 
            (units.erg/units.s, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_luminosity_at_zone", 
            (handler.INDEX, handler.NO_UNIT, units.erg/units.s,), 
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_entropy_at_zone", 
            (handler.INDEX,handler.NO_UNIT,), 
            (units.erg/units.K, handler.ERROR_CODE,)
        )        
        handler.add_method(
            "get_thermal_energy_at_zone", 
            (handler.INDEX,handler.NO_UNIT,), 
            (units.erg/units.g, handler.ERROR_CODE,)
        )        
        handler.add_method(
            "get_id_of_species", 
            (handler.INDEX,handler.NO_UNIT,), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_mass_of_species", 
            (handler.INDEX,handler.NO_UNIT,), 
            (units.amu, handler.ERROR_CODE,)
        )
        handler.add_method(
            "erase_memory", 
            (handler.INDEX,), 
            (handler.ERROR_CODE,),
            public_name = "_erase_memory"
        )
        handler.add_method(
            "new_stellar_model", 
            (units.MSun, units.cm, units.g / units.cm**3, units.K, units.erg / units.s, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,), 
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "finalize_stellar_model", 
            (units.yr,), 
            (handler.INDEX, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_max_age_stop_condition", 
            (), 
            (units.yr, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_max_age_stop_condition", 
            (units.yr, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_min_timestep_stop_condition", 
            (), 
            (units.s, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_min_timestep_stop_condition", 
            (units.s, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_max_iter_stop_condition", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_max_iter_stop_condition", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_mixing_length_ratio", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_mixing_length_ratio", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_semi_convection_efficiency", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_semi_convection_efficiency", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_RGB_wind_scheme", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_RGB_wind_scheme", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_AGB_wind_scheme", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_AGB_wind_scheme", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_reimers_wind_efficiency", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_reimers_wind_efficiency", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_blocker_wind_efficiency", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_blocker_wind_efficiency", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_de_jager_wind_efficiency", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_de_jager_wind_efficiency", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "get_dutch_wind_efficiency", 
            (), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        
    
        handler.add_method(
            "set_dutch_wind_efficiency", 
            (handler.NO_UNIT, ), 
            (handler.ERROR_CODE,)
        )
        
    
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
        self.initialize_code()
        
    def initialize_module_with_current_parameters(self):
        self.initialize_code()
    
    def commit_parameters(self):
        self.parameters.send_not_set_parameters_to_code()
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
        
    def get_mass_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying mass profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_mass_fraction_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_cumulative_mass_profile(self, indices_of_the_stars, number_of_zones = None):
        frac_profile = self.get_mass_profile(indices_of_the_stars, number_of_zones = number_of_zones)
        return frac_profile.cumsum()
    
    def set_mass_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting mass profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_mass_fraction_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
        if hasattr(self, "_erase_memory"):
            self._erase_memory(indices_of_the_stars)
    
    def get_luminosity_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying luminosity profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_luminosity_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def set_luminosity_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting luminosity profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_luminosity_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
        if hasattr(self, "_erase_memory"):
            self._erase_memory(indices_of_the_stars)

    def get_entropy_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying entropy profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_entropy_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_thermal_energy_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying thermal energy profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_thermal_energy_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)

    def get_brunt_vaisala_frequency_squared_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying brunt-vaisala-frequency-squared profiles") 
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_brunt_vaisala_frequency_squared_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_IDs_of_species(self, indices_of_the_stars, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance IDs")
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        return list(self.get_id_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1) 
        ))
    
    def get_masses_of_species(self, indices_of_the_stars, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance mass numbers")
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        return self.get_mass_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1)
        )
    
    def new_particle_from_model(self, internal_structure, current_age=0|units.Myr, key=None):
        if isinstance(internal_structure, dict):
            if "dmass" in internal_structure:
                mass_profile = internal_structure['dmass'][::-1]
            else:
                cumulative_mass_profile = [0.0] | units.MSun
                cumulative_mass_profile.extend(internal_structure['mass'])
                mass_profile = (cumulative_mass_profile[1:] - cumulative_mass_profile[:-1])[::-1]
            self.new_stellar_model(
                mass_profile,
                internal_structure['radius'][::-1],
                internal_structure['rho'][::-1],
                internal_structure['temperature'][::-1],
                internal_structure['luminosity'][::-1],
                internal_structure['X_H'][::-1],
                internal_structure['X_He'][::-1],
                internal_structure['X_C'][::-1],
                internal_structure['X_N'][::-1],
                internal_structure['X_O'][::-1],
                internal_structure['X_Ne'][::-1],
                internal_structure['X_Mg'][::-1],
                internal_structure['X_Si'][::-1],
                internal_structure['X_Fe'][::-1]
            )
        else:
            if hasattr(internal_structure, "dmass"):
                mass_profile = internal_structure.dmass[::-1]
            else:
                cumulative_mass_profile = [0.0] | units.MSun
                cumulative_mass_profile.extend(internal_structure.mass)
                mass_profile = (cumulative_mass_profile[1:] - cumulative_mass_profile[:-1])[::-1]
            self.new_stellar_model(
                mass_profile,
                internal_structure.radius[::-1],
                internal_structure.rho[::-1],
                internal_structure.temperature[::-1],
                internal_structure.luminosity[::-1],
                internal_structure.X_H[::-1],
                internal_structure.X_He[::-1],
                internal_structure.X_C[::-1],
                internal_structure.X_N[::-1],
                internal_structure.X_O[::-1],
                internal_structure.X_Ne[::-1],
                internal_structure.X_Mg[::-1],
                internal_structure.X_Si[::-1],
                internal_structure.X_Fe[::-1]
            )
        tmp_star = datamodel.Particle(key=key)
        tmp_star.age_tag = current_age
        return self.imported_stars.add_particle(tmp_star)

