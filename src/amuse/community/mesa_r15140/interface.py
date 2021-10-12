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
        
    @remote_function   
    def get_maximum_number_of_stars():
        """
        Retrieve the maximum number of stars that can be
        handled by this instance.
        """
        returns (maximum_number_of_stars='i')


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
        
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('filename', dtype='string', direction=function.IN
            , description="The filename of the model to load")
        function.result_type = 'int32'
        return function
        
    @remote_function
    def set_time_step(index_of_the_star='i',time_step='d' | units.yr):
        returns ()
    
    
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
    
    @remote_function   
    def get_mass_loss_rate(index_of_the_star='i'):
        """
        Retrieve the current mass loss rate of the star. (positive for winds, negative for accretion)
        """
        returns (mass_change='d'| units.MSun/units.yr)
    
    @remote_function  
    def get_manual_mass_transfer_rate(index_of_the_star='i'):
        """
        Retrieve the current user-specified mass transfer rate of the star. (negative for winds, positive for accretion)
        """
        returns (mass_change='d'| units.MSun/units.yr)
    
    @remote_function
    def set_manual_mass_transfer_rate(index_of_the_star='i', mass_change='d'| units.MSun/units.yr):
        """
        Set a new user-specified mass transfer rate of the star. (negative for winds, positive for accretion)
        """
        return ()
    
    @remote_function 
    def get_mass_fraction_at_zone(index_of_the_star='i',zone='i'):
        """
        Retrieve the mass fraction at the specified zone/mesh-cell of the star.
        """
        returns (dq_i='d')
    
    
    @remote_function 
    def get_luminosity_at_zone(index_of_the_star='i',zone='i'):
        """
        Retrieve the luminosity at the specified zone/mesh-cell of the star.
        """
        returns (L_i='d')
    
    @remote_function 
    def get_entropy_at_zone(index_of_the_star='i',zone='i'):
        """
        Retrieve the entropy at the specified zone/mesh-cell of the star.
        """
        returns (S_i='d')

    @remote_function 
    def get_thermal_energy_at_zone(index_of_the_star='i',zone='i'):
        """
        Retrieve the thermal energy at the specified zone/mesh-cell of the star.
        """
        returns (E_i='d') 

    @remote_function 
    def get_brunt_vaisala_frequency_squared_at_zone(index_of_the_star='i',zone='i'):
        """
        Retrieve the mass fraction at the specified zone/mesh-cell of the star.
        """
        returns (brunt_N2='d')
    
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
        function.can_handle_array = True
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
        function.can_handle_array = True
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
        function.can_handle_array = True
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
    

    @remote_function 
    def get_control_dble(index_of_the_star='i',name='s'):
        """
        Retrieve the current control option given by name, if it is a double precision number
        """
        returns (value='d')
    
    @remote_function 
    def set_control_dble(index_of_the_star='i',name='s',value='d'):
        """
        Set the current control option given by name, if it is a double precision number
        """
        returns ()

    @remote_function 
    def get_control_logical(index_of_the_star='i',name='s'):
        """
        Retrieve the current control option given by name, if it is a logical flag
        """
        returns (value='b')
    
    @remote_function 
    def set_control_logical(index_of_the_star='i',name='s',value='b'):
        """
        Set the current control option given by name, if it is a logical flag
        """
        returns ()

    @remote_function 
    def get_control_int(index_of_the_star='i',name='s'):
        """
        Retrieve the current control option given by name, if it is an integer
        """
        returns (value='i')
    
    @remote_function 
    def set_control_int(index_of_the_star='i',name='s',value='i'):
        """
        Set the current control option given by name, if it is an integer
        """
        returns ()

    @remote_function 
    def get_control_str(index_of_the_star='i',name='s'):
        """
        Retrieve the current control option given by name, if it is a string
        """
        returns (value='s')
    
    @remote_function 
    def set_control_str(index_of_the_star='i',name='s',value='s'):
        """
        Set the current control option given by name, if it is a string
        """
        returns ()


    @remote_function 
    def get_star_job_dble(index_of_the_star='i',name='s'):
        """
        Retrieve the current star_job option given by name, if it is a double precision number
        """
        returns (value='d')
    
    @remote_function 
    def set_star_job_dble(index_of_the_star='i',name='s',value='d'):
        """
        Set the current star_job option given by name, if it is a double precision number
        """
        returns ()

    @remote_function 
    def get_star_job_logical(index_of_the_star='i',name='s'):
        """
        Retrieve the current star_job option given by name, if it is a logical flag
        """
        returns (value='b')
    
    @remote_function 
    def set_star_job_logical(index_of_the_star='i',name='s',value='b'):
        """
        Set the current star_job option given by name, if it is a logical flag
        """
        returns ()

    @remote_function 
    def get_star_job_int(index_of_the_star='i',name='s'):
        """
        Retrieve the current star_job option given by name, if it is an integer
        """
        returns (value='i')
    
    @remote_function 
    def set_star_job_int(index_of_the_star='i',name='s',value='i'):
        """
        Set the current star_job option given by name, if it is an integer
        """
        returns ()

    @remote_function
    def get_star_job_str(index_of_the_star='i',name='s'):
        """
        Retrieve the current star_job option given by name, if it is a string
        """
        returns (value='s') 
    
    @remote_function
    def set_star_job_str(index_of_the_star='i',name='s',value='s'):
        """
        Set the current star_job option given by name, if it is a string
        """
        returns () 

    @legacy_function   
    def new_stellar_model():
        """
        Define a new star model in the code. 
        Arrays should be in MESA order (surface =1 center=n)
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for par in ['star_mass','dq', 'rho', 'temperature', 
                'XH1','XHe3','XHe4','XC12','XN14','XO16','XNe20','XMg24','XSi28','XS32',
                'XAr36','XCa40','XTi44','XCr48','XFe52','XFe54','XFe56','XCr56','XNi56',
                'prot','neut']:
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


    @remote_function
    def get_profile_at_zone(index_of_the_star='i',zone='i',name='s'):
        """
        Retrieve arbitary profile column at the specified zone/mesh-cell of the star.
        """
        returns (value='d') 

    @remote_function
    def get_history(index_of_the_star='i',name='s'):
        """
        Retrieve arbitary history column of the star.
        """
        returns (value='d')


    @remote_function
    def star_job_update(index_of_the_star='i'):
        """
        After changing options in star_job this function must be called to make the changes to the star
        """
        returns ()


    @remote_function
    def get_nuclear_network(index_of_the_star='i'):
        """
        Retrieve the current nuclear network of the star.
        """
        returns (net_name='s')

    @remote_function
    def set_nuclear_network(index_of_the_star='i',net_name='s'):
        """
        Set the current nuclear network of the star.
        """
        returns ()

    @remote_function   
    def get_age(index_of_the_star='i'):
        """
        Retrieve the current age of the star
        """
        returns (net_name='d')

    @remote_function  
    def set_age(index_of_the_star='i',new_age='d'):
        """
        Set the current age of the star
        """
        returns ()


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

        star_types = ['native_stars','pre_ms_stars',
                    'pre_built_stars','pure_he_stars',
                    'imported_stars'
                    ]

        handler.define_super_set('particles', star_types, 
            index_to_default_set = 0)
        
        handler.define_set('imported_stars', 'index_of_the_star')
        handler.set_new('imported_stars', 'finalize_stellar_model')
        handler.set_delete('imported_stars', 'delete_star')
        
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

        
        for particle_set_name in star_types:
            handler.add_getter(particle_set_name, 'get_radius', names = ('radius',))
            handler.add_getter(particle_set_name, 'get_stellar_type', names = ('stellar_type',))
            handler.add_getter(particle_set_name, 'get_mass', names = ('mass',))
            handler.add_setter(particle_set_name, 'set_mass', names = ('mass',))
            handler.add_getter(particle_set_name, 'get_core_mass', names = ('core_mass',))
            handler.add_getter(particle_set_name, 'get_mass_loss_rate', names = ('wind',))
            handler.add_getter(particle_set_name, 'get_age', names = ('age',))
            handler.add_getter(particle_set_name, 'set_age', names = ('age',))
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
            (units.MSun, handler.NO_UNIT, units.g / units.cm**3, units.K,  
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, 
                handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT,
                handler.NO_UNIT,
                ), 
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_star_age",
            (handler.INDEX,),
            (units.yr, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_star_age",
            (handler.INDEX,units.yr),
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
            "set_age",
            (handler.INDEX,),
            (units.yr,handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_age",
            (handler.INDEX,units.yr),
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

        handler.add_method(
            "start_composition_change",
            (handler.INDEX,handler.NO_UNIT,handler.NO_UNIT,handler.NO_UNIT),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "end_composition_change",
            (handler.INDEX),
            (handler.ERROR_CODE,),
        )

        handler.add_method(
            "set_composition_vector",
            (handler.INDEX,handler.NO_UNIT,handler.NO_UNIT,handler.NO_UNIT),
            (handler.ERROR_CODE,),
        )

        handler.add_method(
            "finalize_stellar_model", 
            (units.yr,), 
            (handler.INDEX, handler.ERROR_CODE,)
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

    def new_particle_from_model(self, internal_structure, current_age=0|units.yr, key=None):
        
        
        if hasattr(internal_structure,'get_profile'):
            mm = internal_structure

            xqs = 10**mm.get_profile('logxq') # log(1-q)
            temperature = mm.get_profile('temperature') | units.K
            rho = mm.get_profile('rho') | units.g / units.cm**3
            
            empty_comp = numpy.zeros(numpy.size(xqs))
            empty_comp[:] = 10**-50 # If you dont have anything in that isotope

            sm = mm.get_history('star_mass')
            star_mass = [sm]*len(xqs) | units.MSun
            species = [mm.get_name_of_species(i) for i in range(1,mm.get_number_of_species()+1)]

            xH1 = mm.get_profile('h1') if 'h1' in species else empty_comp
            xHe3 = mm.get_profile('he3') if 'he3' in species else empty_comp
            xHe4 = mm.get_profile('he4') if 'he4' in species else empty_comp
            xC12 = mm.get_profile('c12') if 'c12' in species else empty_comp
            xN14 = mm.get_profile('n14') if 'n14' in species else empty_comp
            xO16 = mm.get_profile('o16') if 'o16' in species else empty_comp
            xNe20 = mm.get_profile('ne20') if 'ne20' in species else empty_comp
            xMg24 = mm.get_profile('mg24') if 'mg24' in species else empty_comp
            xSi28 = mm.get_profile('si28') if 'si28' in species else empty_comp
            xS32 = mm.get_profile('s32') if 's32' in species else empty_comp
            xAr36 = mm.get_profile('ar36') if 'ar36' in species else empty_comp
            xCa40 = mm.get_profile('ca40') if 'ca40' in species else empty_comp
            xTi44 = mm.get_profile('ti44') if 'ti44' in species else empty_comp
            xCr48 = mm.get_profile('cr48') if 'cr48' in species else empty_comp
            xFe52 = mm.get_profile('fe52') if 'fe52' in species else empty_comp
            xFe54 = mm.get_profile('fe54') if 'fe54' in species else empty_comp
            xFe56 = mm.get_profile('fe56') if 'fe56' in species else empty_comp
            xCr56 = mm.get_profile('cr56') if 'cr56' in species else empty_comp
            xNi56 = mm.get_profile('ni56') if 'ni56' in species else empty_comp

            self.new_stellar_model(
                star_mass,
                xqs[::-1],
                rho[::-1],
                temperature[::-1],
                xH1[::-1],
                xHe3[::-1],
                xHe4[::-1],
                xC12[::-1],
                xN14[::-1],
                xO16[::-1],
                xNe20[::-1],
                xMg24[::-1],
                xSi28[::-1],
                xS32[::-1],
                xAr36[::-1],
                xCa40[::-1],
                xTi44[::-1],
                xCr48[::-1],
                xFe52[::-1],
                xFe54[::-1],
                xFe56[::-1],
                xCr56[::-1],
                xNi56[::-1],
                empty_comp[::-1],
                empty_comp[::-1],
            )


        
        elif isinstance(internal_structure, dict):
            if "dq" in internal_structure:
                mass_profile = internal_structure['dq'][::-1]
            else:
                cumulative_mass_profile = [0.0] | units.MSun
                cumulative_mass_profile.extend(internal_structure['mass'])
                mass_profile = (cumulative_mass_profile[1:] - cumulative_mass_profile[:-1])[::-1]
                mass_profile = mass_profile/cumulative_mass_profile[-1]

            star_mass = numpy.zeros(numpy.size(internal_structure['temperature']))
            star_mass[:] = cumulative_mass_profile[-1]
            empty_comp = numpy.zeros(numpy.size(internal_structure['temperature']))
            empty_comp[:] = 10**-50

            self.new_stellar_model(
                star_mass,
                mass_profile,
                internal_structure['rho'][::-1],
                internal_structure['temperature'][::-1],
                internal_structure['X_H'][::-1],
                empty_comp,
                internal_structure['X_He'][::-1],
                internal_structure['X_C'][::-1],
                internal_structure['X_N'][::-1],
                internal_structure['X_O'][::-1],
                internal_structure['X_Ne'][::-1],
                internal_structure['X_Mg'][::-1],
                internal_structure['X_Si'][::-1],
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                internal_structure['X_Fe'][::-1],
                empty_comp,
                empty_comp,
                empty_comp[::-1],
                empty_comp[::-1],
            )
        else:
            if hasattr(internal_structure, "dq"):
                mass_profile = internal_structure.dq[::-1]
            else:
                cumulative_mass_profile = [0.0] | units.MSun
                cumulative_mass_profile.extend(internal_structure.mass)
                mass_profile = (cumulative_mass_profile[1:] - cumulative_mass_profile[:-1])[::-1]
                mass_profile = mass_profile/cumulative_mass_profile[-1]

            star_mass = numpy.zeros(numpy.size(internal_structure.temperature))
            star_mass[:] = cumulative_mass_profile[-1]
            empty_comp = numpy.zeros(numpy.size(internal_structure.temperature))
            empty_comp[:] = 10**-50


            self.new_stellar_model(
                star_mass,
                mass_profile,
                internal_structure.rho[::-1],
                internal_structure.temperature[::-1],
                internal_structure.X_H[::-1],
                empty_comp,
                internal_structure.X_He[::-1],
                internal_structure.X_C[::-1],
                internal_structure.X_N[::-1],
                internal_structure.X_O[::-1],
                internal_structure.X_Ne[::-1],
                internal_structure.X_Mg[::-1],
                internal_structure.X_Si[::-1],
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                empty_comp,
                internal_structure.X_Fe[::-1],
                empty_comp,
                empty_comp,
                empty_comp[::-1],
                empty_comp[::-1],
            )
        tmp_star = datamodel.Particle(key=key)
        tmp_star.age_tag = current_age
        return self.imported_stars.add_particle(tmp_star)


Mesa = MESA
