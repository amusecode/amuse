from amuse.legacy import *
from amuse.legacy.interface.se import StellarEvolution
from amuse.legacy.support.lit import LiteratureRefs
from amuse.support.interface import CodeInterface
import os

class MESAInterface(LegacyInterface, LiteratureRefs, StellarEvolution): 
    """
        .. [#]  MESA is free software; you can redistribute it and/or modify
                it under the terms of the GNU General Library Public License.
    """
    def __init__(self):
        try:
            LegacyInterface.__init__(self, name_of_the_worker="worker_code")
            LiteratureRefs.__init__(self)
            self.MESA_exists = True
        except Exception:
            print "MESA was not built. Skipping initialization."
            self.MESA_exists = False

    @property
    def default_path_to_inlist(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'mesa_reqs', 'AMUSE_inlist')

    @property
    def default_path_to_MESA_data(self):
        dir = os.path.dirname(__file__)
        return os.path.join(dir, 'src', 'data')

    @legacy_function
    def set_MESA_paths():
        """
        Set the paths to the MESA inlist and data directory.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('inlist_path', dtype='string', direction=function.IN,
            description = "Path to the inlist file.")
        function.addParameter('data_path', dtype='string', direction=function.IN,
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

class MESA(CodeInterface):
    
    def __init__(self):
        CodeInterface.__init__(self, MESAInterface())
        if self.MESA_exists:
            self.set_MESA_paths(self.default_path_to_inlist, 
                self.default_path_to_MESA_data)
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
        
        object.add_method('particles', 'evolve', 'evolve_one_step')
    
    def define_errorcodes(self, object):
        object.add_errorcode(-1, 'Something went wrong...')
        object.add_errorcode(-11, 'Evolve terminated: Unspecified stop condition reached.')
        object.add_errorcode(-12, 'Evolve terminated: Maximum age reached.')
        object.add_errorcode(-13, 'Evolve terminated: Maximum number of iterations reached.')
    
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
        
    
    def initialize_module_with_default_parameters(self):
        self.parameters.set_defaults()
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
                
        
