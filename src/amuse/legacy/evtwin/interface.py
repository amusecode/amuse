from amuse.legacy import *
from amuse.legacy.interface.se import StellarEvolution
from amuse.legacy.support.lit import LiteratureRefs
from amuse.support.interface import CodeInterface

import os

class EVtwinInterface(LegacyInterface, LiteratureRefs, StellarEvolution):
    """
    Need to have docs
    """
    use_modules = ['twin_library_v2']
    
    def __init__(self):
        LegacyInterface.__init__(self, name_of_the_worker="worker_code")
        LiteratureRefs.__init__(self)
         
    @property
    def default_path_to_ev_database(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'src')
        
    @legacy_function   
    def new_particle():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('mass', dtype='float64', direction=function.IN
            , description="The initial mass of the star")
        
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            New star could not be created.
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
    def set_maximum_number_of_stars():
        """
        Update the maximum number of stars that can be
        handled by this instance. Need to set this number
        before calling :method:`initialize_code`. Cannot be
        changed once initialize_code has been called.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('maximum_number_of_stars', dtype='int32', direction=function.IN,
            description = "The new value of the maximum number of stars.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            The code cannot update the maximum number of stars
        """
        return function
        
    
    @legacy_function
    def set_ev_path():
        """
        Update the path to the EVtwin database.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='string', direction=function.IN,
            description = "Name of the the directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            Directory does not exist
        """
        return function
        
    @legacy_function
    def set_init_dat_name():
        """
        Update name of the init.dat file
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='string', direction=function.IN,
            description = "File in the evpath directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
        """
        return function
        
    @legacy_function
    def set_init_run_name():
        """
        Update name of the init.run file
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('path', dtype='string', direction=function.IN,
            description = "File in the evpath directory")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was set
        -1 - ERROR
            File does not exist
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
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
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
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
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
    def get_convective_overshoot_parameter():
        """
        Retrieve the current value of the convective overshoot parameter.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('convective_overshoot_parameter', dtype='float64', direction=function.OUT
            , description="The current value of the convective overshoot parameter.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_convective_overshoot_parameter():
        """
        Set the value of the convective overshoot parameter.
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('convective_overshoot_parameter', dtype='float64', direction=function.IN
            , description="The new value of the convective overshoot parameter.")
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
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
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
        after Langer, Sugimoto & Fricke 1983 (A&A).
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
        after Langer, Sugimoto & Fricke 1983 (A&A).
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
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
    def get_thermohaline_mixing_parameter():
        """
        Retrieve the current value of the thermohaline mixing parameter,
        probably only important for binaries and collision remnants.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('thermohaline_mixing_parameter', dtype='float64', direction=function.OUT
            , description="The current value of the thermohaline mixing parameter.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_thermohaline_mixing_parameter():
        """
        Set the value of the thermohaline mixing parameter,
        probably only important for binaries and collision remnants.
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('thermohaline_mixing_parameter', dtype='float64', direction=function.IN
            , description="The new value of the thermohaline mixing parameter.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_number_of_ionization_elements():
        """
        Retrieve the current number of elements used for ionization of this instance.
        With the default value (2), only the ionization of H and He are taken into account
        in the EoS. For values 3, 4, 5, 6, 7, 8, 9:
        the elements:          C, N, O, Ne,Mg,Si,Fe are also included. Don't try 9.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_of_ionization_elements', dtype='int32', direction=function.OUT
            , description="The current number of elements used for ionization in EoS solver of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_number_of_ionization_elements():
        """
        Set the new number of elements used for ionization of this instance.
        With the default value (2), only the ionization of H and He are taken into account
        in the EoS. For values 3, 4, 5, 6, 7, 8, 9:
        the elements:          C, N, O, Ne,Mg,Si,Fe are also included. Don't try 9.
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_of_ionization_elements', dtype='int32', direction=function.IN
            , description="The new number of elements used for ionization in EoS solver of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_AGB_wind_setting():
        """
        Retrieve the current AGB wind setting of this instance (1 or 2).
        (1) means use Wachter et al. (AGB) mass loss prescription.
        (2) means use Vasiliadis & Wood (AGB) mass loss rate.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('AGB_wind_setting', dtype='int32', direction=function.OUT
            , description="The current AGB wind setting of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_AGB_wind_setting():
        """
        Set the new AGB wind setting of this instance (1 or 2).
        (1) means use Wachter et al. (AGB) mass loss prescription.
        (2) means use Vasiliadis & Wood (AGB) mass loss rate.
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('AGB_wind_setting', dtype='int32', direction=function.IN
            , description="The new AGB wind setting of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function
        
    @legacy_function
    def get_RGB_wind_setting():
        """
        Retrieve the current RGB wind setting of this instance.
        (positive) means use Schroeder & Cuntz mass loss prescription.
        (negative) means use Reimers mass loss rate. (0) means none.
        The absolute value is used as efficiency factor.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('RGB_wind_setting', dtype='float64', direction=function.OUT
            , description="The current RGB wind setting of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
    
    @legacy_function
    def set_RGB_wind_setting():
        """
        Set the new RGB wind setting of this instance.
        (positive) means use Schroeder & Cuntz mass loss prescription.
        (negative) means use Reimers mass loss rate. (0) means none.
        The absolute value is used as efficiency factor.
        This needs to be set after calling :method:`initialize_code`. It will 
        be overridden by initialize_code otherwise.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('RGB_wind_setting', dtype='float64', direction=function.IN
            , description="The new RGB wind setting of this instance.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            The code could not set the value.
        """
        return function

    @legacy_function
    def get_time_step():
        """
        Retrieve the current time step (yr) to be taken for the evolution of this star.
        Note that the stellar evolution code might change the value during an 
        evolve_model call, if it fails to converge using the current value.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the stellar type of")
        function.addParameter('time_step', dtype='float64', direction=function.OUT
            , description="The current time step (yr) to be taken for the evolution of this star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function

    @legacy_function   
    def new_spinning_particle():
        """
        Define a new star in the code. The star will start with the given mass and given spin.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('mass', dtype='float64', direction=function.IN
            , description="The initial mass of the star")
        function.addParameter('spin', dtype='float64', direction=function.IN
            , description="The initial spin of the star")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            New star could not be created.
        """
        return function
        
    @legacy_function
    def get_spin():
        """
        Retrieve the current spin period (in days) of this star.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the stellar type of")
        function.addParameter('spin', dtype='float64', direction=function.OUT
            , description="The current spin period (in days) of this star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code could not retrieve the value.
        """
        return function
        
class EVtwin(CodeInterface):
    
    def __init__(self):
        CodeInterface.__init__(self, EVtwinInterface())
        
    
    def define_parameters(self, object):
              
        object.add_method_parameter(
            "get_maximum_number_of_stars",
            "set_maximum_number_of_stars",
            "maximum_number_of_stars", 
            "Maximum number of stars that can be allocated", 
            units.none, 
            10 | units.none
        )
        
        object.add_method_parameter(
            "get_metallicity",
            "set_metallicity",
            "metallicity", 
            "Metallicity of all stats", 
            units.none, 
            0.02 | units.none
        )
        
        object.add_method_parameter(
            None,
            "set_ev_path",
            "path_to_data", 
            "Path to the data directory", 
            units.string, 
            None
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
            "get_number_of_ionization_elements",
            "set_number_of_ionization_elements",
            "number_of_ionization_elements", 
            "The number of elements used for ionization in EoS solver of this instance.",
            units.none, 
            2 | units.none
        )
        
        object.add_method_parameter(
            "get_convective_overshoot_parameter",
            "set_convective_overshoot_parameter",
            "convective_overshoot_parameter", 
            "The convective overshoot parameter.",
            units.none, 
            0.12 | units.none
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
            "The efficiency of semi-convection, after Langer, Sugimoto & Fricke 1983 (A&A).",
            units.none, 
            0.04 | units.none
        )
        
        object.add_method_parameter(
            "get_thermohaline_mixing_parameter",
            "set_thermohaline_mixing_parameter",
            "thermohaline_mixing_parameter", 
            "The thermohaline mixing parameter, probably only important for binaries and collision remnants.",
            units.none, 
            1.0 | units.none
        )
        
        object.add_method_parameter(
            "get_AGB_wind_setting",
            "set_AGB_wind_setting",
            "AGB_wind_setting", 
            "The AGB wind setting: (1, 2) for (Wachter&al, Vasiliadis&Wood) mass loss.",
            units.none, 
            1 | units.none
        )
        
        object.add_method_parameter(
            "get_RGB_wind_setting",
            "set_RGB_wind_setting",
            "RGB_wind_setting", 
            "The RGB wind setting: (positive, negative, 0) for (Schroeder&Cuntz, Reimers, none) mass loss.",
            units.none, 
            1.0 | units.none
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
        object.add_getter('particles', 'get_spin', names = ('spin',))
        object.add_getter('particles', 'get_luminosity',names = ('luminosity',))
        
        object.add_method('particles', 'evolve', 'evolve_one_step')
    
    def define_errorcodes(self, object):
        object.add_errorcode(-2, 'BEGINN -- requested mesh too large')
        object.add_errorcode(-1, 'STAR12 -- no timesteps required')
        object.add_errorcode(2, 'BACKUP -- tstep reduced below limit; quit')
        object.add_errorcode(3, 'NEXTDT -- *2 evolving beyond last *1 model')
        object.add_errorcode(4, 'PRINTB -- *1 rstar exceeds rlobe by limit')
        object.add_errorcode(5, 'PRINTB -- age greater than limit')
        object.add_errorcode(6, 'PRINTB -- C-burning exceeds limit')
        object.add_errorcode(7, 'PRINTB -- *2 radius exceeds rlobe by limit')
        object.add_errorcode(8, 'PRINTB -- close to He flash')
        object.add_errorcode(9, 'PRINTB -- massive (>1.2 msun) deg. C/O core')
        object.add_errorcode(10, 'PRINTB -- |M1dot| exceeds limit')
        object.add_errorcode(11, 'NEXTDT -- impermissible FDT for *2')
        object.add_errorcode(14, 'PRINTB -- funny compos. distribution')
        object.add_errorcode(15, 'STAR12 -- terminated by hand')
        object.add_errorcode(16, 'MAIN   -- ZAHB didnt converge')
        object.add_errorcode(17, 'BEGINN -- Nucleosynthesis did not converge')
        object.add_errorcode(51, 'PRINTB -- end of MS (core H below limit)')
        object.add_errorcode(52, 'PRINTB -- Radius exceeds limit')
        object.add_errorcode(53, 'PRINTB -- Convergence to target model reached minimum')
        object.add_errorcode(12, 'BACKUP -- tstep reduced below limit; quit -- core H non-zero')
        object.add_errorcode(22, 'BACKUP -- tstep reduced below limit; quit -- core He non-zero')
        object.add_errorcode(32, 'BACKUP -- tstep reduced below limit; quit -- core C non-zero')
        
    
    def define_methods(self, object):
            
        object.add_method(
            'evolve', 
            (object.INDEX,), 
            (object.ERROR_CODE,), 
        )
        object.add_method(
            "new_particle", 
            (units.MSun), #, units.day, units.RSun),
            (object.NO_UNIT, object.ERROR_CODE)
        )
        object.add_method(
            "delete_star", 
            (object.NO_UNIT,), 
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_mass", 
            (object.NO_UNIT,), 
            (units.MSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_radius", 
            (object.INDEX,), 
            (units.RSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_stellar_type", 
            (object.NO_UNIT,), 
            (units.stellar_type, object.ERROR_CODE,)
        )
        object.add_method(
            "get_age", 
            (object.NO_UNIT,), 
            (units.yr, object.ERROR_CODE,)
        )
        object.add_method(
            "get_luminosity", 
            (object.NO_UNIT,), 
            (units.LSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_time_step", 
            (object.NO_UNIT,), 
            (units.yr, object.ERROR_CODE,)
        )
        object.add_method(
            "get_spin", 
            (object.NO_UNIT,), 
            (units.day, object.ERROR_CODE,)
        )
        
    
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
        self.set_ev_path(self.default_path_to_ev_database)
        self.initialize_code()
        
    def initialize_module_with_current_parameters(self):
        self.initialize_code()
        
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def evolve_model(self, end_time = None):
        if end_time is None:
            result = self.particles.evolve_one_step()
            return result
                   
        for particle in self.particles:
            while particle.age < end_time:
                particle.evolve_one_step()
                
        
