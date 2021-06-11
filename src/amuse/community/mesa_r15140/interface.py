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
    http://mesa.sourceforge.net/) version 15140, aims to provide state-of-the-art, robust, 
    and efficient open source modules, usable singly or in combination for a 
    wide range of applications in stellar astrophysics. The AMUSE interface to 
    MESA can create and evolve stars using the MESA/STAR module. 
    All metallicities are supported, even the 
    interesting case of Z=0. The supported stellar mass range is from 
    about 1M_jupiter to >100 Msun.
    
    References:
        .. [#] Paxton, Bildsten, Dotter, Herwig, Lesaffre & Timmes 2011, ApJS, arXiv:1009.1622 [2011ApJS..192....3P]
        .. [#] Paxton, Cantiello, Arras, Bildsten, Brown, Dotter, Mankovich, Montgomery, Stello, Timmes, Townsend, 2013, ApJS, arXiv:1301.0319, [2013ApJS..208....4P]
        .. [#] Paxton, Marchant, Schwab, Bauer, Bildsten, Cantiello, Dessart, Farmer, Hu, Langer, Townsend, Townsley, Timmes, 2015, ApJS, arXiv:1506.03146, [2015ApJS..220...15P]
        .. [#] Paxton, Schwab, Bauer, Bildsten, Blinnikov, Duffell, Farmer, Goldberg, Marchant, Sorokina, Thoul, Townsend, Timmes, 2018, arXiv:1710.08424, [2018ApJS..234...34P] 
        .. [#] Paxton, Smolec, Schwab, Gautschy, Bildsten, Cantiello, Dotter, Farmer, Goldberg, Jermyn, Kanbur, Marchant, Thoul, Townsend, Wolf, Zhang, Timmes, [2019ApJS..243...10P]
        .. [#] http://mesa.sourceforge.net/
        .. [#] https://docs.mesastar.org/en/latest/reference.html
    """

    use_modules = ['amuse_mesa']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="mesa_worker", **options)
        LiteratureReferencesMixIn.__init__(self)
        CodeWithDataDirectories.__init__(self)
    
    @property
    def default_path_to_inlist(self):
        return ''

    @option(type="string", sections=('data'))
    def default_path_to_MESA(self):
        return os.path.join(self.amuse_root_directory, 'src', 'amuse', 'community', 'mesa_r15140', 'src', 'mesa-r15140')
    

    @legacy_function
    def set_MESA_paths():
        """
        Set the paths to the MESA inlist and data directories.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('inlist_path', dtype='string', direction=function.IN,
            description = "Path to the inlist file.")
        function.addParameter('mesa_dir', dtype='string', direction=function.IN,
            description = "Path to the MESA directory.")
        function.addParameter('local_data_path', dtype='string', direction=function.IN,
            description = "Path to the data directory.")
        function.addParameter('gyre_in_filename', dtype='string', direction=function.IN,
            description = "Path to the gyre.in file.")
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
    def new_pure_he_particle():
        """
        Define a new pure He star in the code. The star will start with the given mass.
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
    def new_zams_particle():
        """
        Define a new ZAMS model with Z=0.02. The star will start with the given mass.
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
    def load_model():
        """
        Load a pre-built MESA model (.mod file)
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('filename', dtype='string', direction=function.IN
            , description="The filename of the model to load")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('time_step', dtype='float64', direction=function.IN
            , description="The next timestep for the star in years.")
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
        Retrieve the net_id of the chemical abundance variable of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='string', direction=function.IN
            , description="The name of the isotope to get the name id of")
        function.addParameter('species_id', dtype='int32', direction=function.OUT
            , description="The net_id of the chemical abundance variable of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            The isotope was not found.
        """
        return function
    
    @legacy_function
    def get_mass_of_species():
        """
        Retrieve the mass number of the species.
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
        -2 - ERROR
            The species was not found
        """
        return function
    
    @legacy_function
    def get_mass_fraction_of_species_at_zone():
        """
        Retrieve the mass number of the chemical abundance variable of the star at zone.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to get the mass number of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone of the star to get the mass number of")
        function.addParameter('species_mass', dtype='float64', direction=function.OUT
            , description="The mass number of the chemical abundance variable of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            The zone was not found
        -3 - ERROR
            The species was not found
        """
        return function

    @legacy_function
    def get_name_of_species():
        """
        Retrieve the name of the species given by the species id
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to get the mass number of")
        function.addParameter('species_name', dtype='string', direction=function.OUT
            , description="The name of the species.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            The species was not found
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
    def get_control_dble():
        """
        Retrieve the current control option given by name, if it is a double precision number
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_control_dble():
        """
        Set the current control option given by name, if it is a double precision number
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='float64', direction=function.IN
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_control_logical():
        """
        Retrieve the current control option given by name, if it is a logical flag
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='bool', direction=function.OUT
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_control_logical():
        """
        Set the current control option given by name, if it is a logical flag
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='bool', direction=function.IN
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_control_int():
        """
        Retrieve the current control option given by name, if it is an integer
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='int32', direction=function.OUT
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_control_int():
        """
        Set the current control option given by name, if it is an integer
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='int32', direction=function.IN
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_control_str():
        """
        Retrieve the current control option given by name, if it is a string
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='string', direction=function.OUT
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_control_str():
        """
        Set the current control option given by name, if it is a string
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the control to get")
        function.addParameter('value', dtype='string', direction=function.IN
            , description="The value of the control option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function


    @legacy_function
    def get_star_job_dble():
        """
        Retrieve the current star_job option given by name, if it is a double precision number
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_star_job_dble():
        """
        Set the current star_job option given by name, if it is a double precision number
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='float64', direction=function.IN
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_star_job_logical():
        """
        Retrieve the current star_job option given by name, if it is a logical flag
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='bool', direction=function.OUT
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_star_job_logical():
        """
        Set the current star_job option given by name, if it is a logical flag
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='bool', direction=function.IN
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_star_job_int():
        """
        Retrieve the current star_job option given by name, if it is an integer
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='int32', direction=function.OUT
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_star_job_int():
        """
        Set the current star_job option given by name, if it is an integer
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='int32', direction=function.IN
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_star_job_str():
        """
        Retrieve the current star_job option given by name, if it is a string
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='string', direction=function.OUT
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_star_job_str():
        """
        Set the current star_job option given by name, if it is a string
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', 
            direction=function.IN, description = "The index of the star. ")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The name of the star_job to get")
        function.addParameter('value', dtype='string', direction=function.IN
            , description="The value of the star_job option")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function




    @legacy_function   
    def new_stellar_model():
        """
        Define a new star model in the code. 
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
    def get_profile_at_zone():
        """
        Retrieve arbitary profile column at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to get the value of")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The profile column name to get the value of")            
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The profile value at the specified zone/mesh-cell of the star.")
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
    def get_history():
        """
        Retrieve arbitary history column of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('name', dtype='string', direction=function.IN
            , description="The history column name to get the value of")            
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The history value of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        -2 - ERROR
            A coloumn with the given name was not found.
        """
        return function  


    @legacy_function
    def star_job_update():
        """
        After changing options in star_job this function must be called to make the changes to the star
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to update")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The changes were appplied.
        -1 - ERROR
            The changes were not appplied.
        """
        return function  


    @legacy_function
    def get_nuclear_network():
        """
        Retrieve the current nuclear network of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('net_name', dtype='string', direction=function.OUT
            , description="The current nuclear network of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function

    @legacy_function
    def set_nuclear_network():
        """
        Set the current nuclear network of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('net_name', dtype='string', direction=function.IN
            , description="The new nuclear network of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was retrieved.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function


class MESA(StellarEvolution, InternalStellarStructure):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, MESAInterface(**options), **options)
        
        output_dir = self.get_output_directory()

        if 'inlist' in options:
            inlist = options['inlist']
            if not os.path.exists(inlist):
                raise ValueError('Named inlist does not exist, maybe in a different folder?')
        else:
            inlist = self.default_path_to_inlist


        if 'gyre_in' in options:
            gyre_in = options['gyre_in']
        else:
            gyre_in  = ''

        self.set_MESA_paths(
            inlist, 
            self.default_path_to_MESA, 
            output_dir,
            gyre_in
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
        
        
    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['native_stars','pre_ms_stars','pre_built_stars','pure_he_stars'], 
            index_to_default_set = 0)
        
        # handler.define_set('imported_stars', 'index_of_the_star')
        # handler.set_new('imported_stars', 'finalize_stellar_model')
        # handler.set_delete('imported_stars', 'delete_star')
        
        handler.define_set('native_stars', 'index_of_the_star')
        handler.set_new('native_stars', 'new_zams_particle')
        handler.set_delete('native_stars', 'delete_star')
        
        handler.define_set('pre_ms_stars', 'index_of_the_star')
        handler.set_new('pre_ms_stars', 'new_pre_ms_particle')
        handler.set_delete('pre_ms_stars', 'delete_star')

        handler.define_set('pre_built_stars', 'index_of_the_star')
        handler.set_new('pre_built_stars', 'load_model')
        handler.set_delete('pre_built_stars', 'delete_star')

        handler.define_set('pure_he_stars', 'index_of_the_star')
        handler.set_new('pure_he_stars', 'new_pure_he_particle')
        handler.set_delete('pure_he_stars', 'delete_star')

        
        for particle_set_name in ['native_stars','pre_ms_stars','pre_built_stars','pure_he_stars']:
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

            handler.add_method(particle_set_name, 'get_control_dble')
            handler.add_method(particle_set_name, 'set_control_dble')
            handler.add_method(particle_set_name, 'get_control_int')
            handler.add_method(particle_set_name, 'set_control_int')
            handler.add_method(particle_set_name, 'get_control_str')
            handler.add_method(particle_set_name, 'set_control_str')
            handler.add_method(particle_set_name, 'get_control_logical')
            handler.add_method(particle_set_name, 'set_control_logical')

            handler.add_method(particle_set_name, 'get_star_job_dble')
            handler.add_method(particle_set_name, 'set_star_job_dble')
            handler.add_method(particle_set_name, 'get_star_job_int')
            handler.add_method(particle_set_name, 'set_star_job_int')
            handler.add_method(particle_set_name, 'get_star_job_str')
            handler.add_method(particle_set_name, 'set_star_job_str')
            handler.add_method(particle_set_name, 'get_star_job_logical')
            handler.add_method(particle_set_name, 'set_star_job_logical')

            handler.add_method(particle_set_name, 'star_job_update')

            handler.add_method(particle_set_name, 'get_nuclear_network')
            handler.add_method(particle_set_name, 'set_nuclear_network')
                        
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
            handler.add_method(particle_set_name, 'get_profile')
            handler.add_method(particle_set_name, 'get_history') 

            handler.add_method(particle_set_name, 'get_name_of_species')
            handler.add_method(particle_set_name, 'get_mass_of_species')
            handler.add_method(particle_set_name, 'get_masses_of_species')
            handler.add_method(particle_set_name, 'get_mass_fraction_of_species_at_zone')
            handler.add_method(particle_set_name, 'get_id_of_species')
            
            handler.add_method(particle_set_name, 'get_IDs_of_species')


            
    def define_state(self, handler):
        StellarEvolution.define_state(self, handler)
        handler.add_method('EDIT', 'new_pre_ms_particle')
        handler.add_method('UPDATE', 'new_pre_ms_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_pre_ms_particle', False)
        # handler.add_method('EDIT', 'finalize_stellar_model')
        # handler.add_method('UPDATE', 'finalize_stellar_model')
        # handler.add_transition('RUN', 'UPDATE', 'finalize_stellar_model', False)
    
    def define_errorcodes(self, handler):
        InternalStellarStructure.define_errorcodes(self, handler)
        handler.add_errorcode(-1, 'Something went wrong...')
        handler.add_errorcode(-4, 'Not implemented.')
        handler.add_errorcode(-11, 'Evolve terminated: Unspecified stop condition reached.')
        handler.add_errorcode(-12, 'Evolve terminated: Maximum age reached.')
        handler.add_errorcode(-13, 'Evolve terminated: Maximum number of iterations reached.')
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
            "new_zams_particle",
            (units.MSun),
            (handler.INDEX, handler.ERROR_CODE)
        )
        handler.add_method(
            "new_pure_he_particle",
            (units.MSun),
            (handler.INDEX, handler.ERROR_CODE)
        )
        handler.add_method(
            "load_model",
            (handler.NO_UNIT),
            (handler.INDEX, handler.ERROR_CODE)
        )
        handler.add_method(
            "set_time_step", 
            (handler.INDEX, units.s), 
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_time_step", 
            (handler.INDEX,), 
            (units.s, handler.ERROR_CODE,)
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
            "get_profile_at_zone", 
            (handler.INDEX,handler.NO_UNIT, handler.NO_UNIT), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )    
        handler.add_method(
            "get_history", 
            (handler.INDEX, handler.NO_UNIT), 
            (handler.NO_UNIT, handler.ERROR_CODE,)
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
            "new_stellar_model", 
            (units.MSun, units.cm, units.g / units.cm**3, units.K, units.erg / units.s, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,), 
            (handler.ERROR_CODE,)
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
            "get_mass_fraction_of_species_at_zone",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.NO_UNIT,handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_name_of_species",
            (handler.INDEX, handler.NO_UNIT,),
            (handler.NO_UNIT,handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_control_dble",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_control_dble",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_control_int",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_control_int",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_control_str",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_control_str",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_control_logical",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_control_logical",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "get_star_job_dble",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_star_job_dble",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_star_job_int",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_star_job_int",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_star_job_str",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_star_job_str",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_star_job_logical",
            (handler.INDEX,handler.NO_UNIT),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_star_job_logical",
            (handler.INDEX, handler.NO_UNIT, handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "star_job_update",
            (handler.INDEX,),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_nuclear_network",
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_nuclear_network",
            (handler.INDEX,handler.NO_UNIT),
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
        return self.get_mass_fraction_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none)
    
    def get_cumulative_mass_profile(self, indices_of_the_stars, number_of_zones = None):
        frac_profile = self.get_mass_profile(indices_of_the_stars, number_of_zones = number_of_zones)
        return frac_profile.cumsum()
    
    def set_mass_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting mass profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_mass_fraction_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none, values)
    
    def get_luminosity_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying luminosity profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_luminosity_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none)
    
    def set_luminosity_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting luminosity profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_luminosity_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none, values)

    def get_entropy_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying entropy profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_entropy_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none)
    
    def get_thermal_energy_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying thermal energy profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_thermal_energy_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none)

    def get_brunt_vaisala_frequency_squared_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying brunt-vaisala-frequency-squared profiles") 
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_brunt_vaisala_frequency_squared_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none)
    
    def get_IDs_of_species(self, indices_of_the_stars, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance IDs")
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        return list(range(1,number_of_species+1))

    def get_profile(self, indices_of_the_stars, name, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_profile_at_zone([indices_of_the_stars]*number_of_zones, list(range(number_of_zones)) | units.none, [name]*number_of_zones)
        
    def get_masses_of_species(self, indices_of_the_stars, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance mass numbers")
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        return self.get_mass_of_species(
            [indices_of_the_stars]*number_of_species, 
            list(range(1,number_of_species+1))
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



Mesa = MESA
