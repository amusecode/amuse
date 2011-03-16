"""
Stellar Dynamics Interface Defintion
"""

import numpy
from amuse.support.units import units
from amuse.support import exceptions
from amuse.support.codes.core import legacy_function, LegacyFunctionSpecification

class StellarEvolution(object): 


    @legacy_function   
    def delete_star():
        """
        Remove the star with the given index from the code.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to remove")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The star has been deleted
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
       
        
    
    @legacy_function   
    def get_luminosity():
        """
        Retrieve the current luminosity of the star.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('luminosity', dtype='float64', direction=function.OUT
            , description="The current luminosit of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
    @legacy_function   
    def get_mass():
        """
        Retrieve the current mass of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('mass', dtype='float64', direction=function.OUT
            , description="The current mass of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
    
    @legacy_function   
    def get_radius():
        """
        Retrieve the current radius of the star.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('radius', dtype='float64', direction=function.OUT
            , description="The current radius of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """

        return function
        
    @legacy_function   
    def get_temperature():
        """
        Retrieve the current temperature of the star.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('temperature', dtype='float64', direction=function.OUT
            , description="The current temperature of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
    @legacy_function   
    def get_age():
        """
        Retrieve the current age of the star.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('age', dtype='float64', direction=function.OUT
            , description="The current age of the star. ")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
    #@legacy_function   
    #def get_type():
    #    """
    #    Retrieve the type of the star. The meaning of the stellar type is defined
    #    by the code. (Difference between stellar type and type must be explained)
    #    """
    #   function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
    #        , description="The index of the star to get the type of")
    #    function.addParameter('type', dtype='int32', direction=function.OUT
    #        , description="The type. ")
    #    function.result_type = 'i'
    #    function.result_doc = """
    #    0 - OK
    #        The value has been set.
    #    -1 - ERROR
    #        A star with the given index was not found.
    #    """
    #return function
        
    @legacy_function   
    def get_stellar_type():
        """
        Retrieve the stellar type of the star. The meaning of the stellar type is defined
        by the code.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the stellar type of")
        function.addParameter('stellar_type', dtype='int32', direction=function.OUT
            , description="The stellar type. ")
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
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
        #function.addParameter('age_tag', dtype='float64', direction=function.IN
        #    , description="Starting age of the star *to be specified exactly*")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            New star could not be created.
        """
        return function
        
    @legacy_function
    def get_number_of_particles():
        """
        Retrieve the total number of particles define  d in the code
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT,
            description = "Count of the particles in the code")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            Count could be determined
         -1 - ERROR
            Unable to determine the count
        """
        return function 
        
    @legacy_function   
    def evolve():
        """
        Evolve the star with the given index one step. The code determines how far in time the star will
        be evolved after this function is finished. See the ``get_age`` function for retrieving
        the current age of the star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star, as returned by the new_zams_star function")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The star was evolved.
        -1 - ERROR
            The evolution could not complete, solution did not converge.
        """
        return function
        
        
    @legacy_function
    def initialize_code():
        """
        Let the code perform initialization actions after all parameters have been set.
        """
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Code is initialized
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention 
        """
        return function  
        
    
    @legacy_function
    def get_metallicity():
        """
        Retrieve the current value of the metallicity. The range of metallicities
        one can simulate will be dependent on the code.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('metallicity', dtype='float64', direction=function.OUT,
            description = "The current value of the metallicity")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the metallicity was retrieved
        -1 - ERROR
            The code does not have support for retrieving the metallicity
        """
        return function
        
    
    @legacy_function
    def set_metallicity():
        """
        Update the value of the metallicity.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('metallicity', dtype='float64', direction=function.IN,
            description = "The new value of the metallicity.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the metallicity was set
        -1 - ERROR
            The code does not have support for updating the metallicity
        """
        return function
        
    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions after all particles have been created. 
        Called before the first evolve call and after the last new_particle call. Not every
        code needs this functionality. And it may be possible to create stars after this
        function has been called
        """
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention 
        """
        return function  
    


class InternalStellarStructureInterface(object): 

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
    def set_temperature_at_zone():
        """
        Set the temperature at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to set the value of")
        function.addParameter('T_i', dtype='float64', direction=function.IN
            , description="The temperature at the specified zone/mesh-cell of the star.")
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
    def set_density_at_zone():
        """
        Set the density at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to set the value of")
        function.addParameter('rho_i', dtype='float64', direction=function.IN
            , description="The density at the specified zone/mesh-cell of the star.")
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
    def set_mass():
        """
        Set the current mass of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('mass', dtype='float64', direction=function.IN
            , description="The current mass of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
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
    def set_radius_at_zone():
        """
        Set the radius at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to set the value of")
        function.addParameter('R_i', dtype='float64', direction=function.IN
            , description="The radius at the specified zone/mesh-cell of the star.")
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
    def set_mass_fraction_of_species_at_zone():
        """
        Set the fractional chemical abundance variable at the specified zone/mesh-cell of the star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('species', dtype='int32', direction=function.IN
            , description="The species of the star to set the value of")
        function.addParameter('zone', dtype='int32', direction=function.IN
            , description="The zone/mesh-cell of the star to set the value of")
        function.addParameter('Xj_i', dtype='float64', direction=function.IN
            , description="The fractional chemical abundance variable at the specified zone/mesh-cell of the star.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value was set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
    


class InternalStellarStructure(object): 

    def define_particle_sets(self, object, set_name = 'particles'):
        object.add_method(set_name, 'get_number_of_zones')
        object.add_method(set_name, 'get_density_profile')
        object.add_method(set_name, 'set_density_profile')
        object.add_method(set_name, 'get_radius_profile')
        object.add_method(set_name, 'set_radius_profile')
        object.add_method(set_name, 'get_temperature_profile')
        object.add_method(set_name, 'set_temperature_profile')
        object.add_method(set_name, 'get_mu_profile')
        object.add_method(set_name, 'get_number_of_species')
        object.add_method(set_name, 'get_names_of_species')
        object.add_method(set_name, 'get_chemical_abundance_profiles')
        object.add_method(set_name, 'set_chemical_abundance_profiles')
    
    def define_errorcodes(self, object):
        object.add_errorcode(-21, 'Specified particle does not exist.')
        object.add_errorcode(-22, 'Specified zone is undefined for this particle.')
        object.add_errorcode(-23, 'Specified chemical species is undefined for this particle.')
    
    def define_methods(self, object):
        object.add_method(
            "set_mass",
            (object.INDEX, units.MSun,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_number_of_zones", 
            (object.INDEX,), 
            (units.none, object.ERROR_CODE,)
        )
        object.add_method(
            "get_temperature_at_zone", 
            (object.INDEX,units.none,), 
            (units.K, object.ERROR_CODE,)
        )
        object.add_method(
            "set_temperature_at_zone", 
            (object.INDEX, units.none, units.K,), 
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_density_at_zone", 
            (object.INDEX,units.none,), 
            (units.g/units.cm**3, object.ERROR_CODE,)
        )
        object.add_method(
            "set_density_at_zone", 
            (object.INDEX, units.none, units.g/units.cm**3,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_radius_at_zone", 
            (object.INDEX,units.none,), 
            (units.cm, object.ERROR_CODE,)
        )
        object.add_method(
            "set_radius_at_zone", 
            (object.INDEX, units.none, units.cm,), 
            (object.ERROR_CODE,)
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
            "get_mass_fraction_of_species_at_zone", 
            (object.INDEX,units.none,units.none,), 
            (units.none, object.ERROR_CODE,)
        )
        object.add_method(
            "set_mass_fraction_of_species_at_zone", 
            (object.INDEX, units.none, units.none, units.none,), 
            (object.ERROR_CODE,)
        )
    
    def _check_supplied_values(self, number_of_values, expected_number, type_string = "mesh zones"):
        if number_of_values != expected_number:
            raise exceptions.CodeException(("The length of the supplied vector ({0}) does not match the number of "
                +type_string+" of the star ({1}).").format(number_of_values, expected_number))
    
    def _check_number_of_indices(self, indices_of_the_stars, action_string = "Querying/setting profiles"):
        if hasattr(indices_of_the_stars, '__iter__'):
            if len(indices_of_the_stars) > 1:
                raise exceptions.CodeException(action_string+" of more than one particle at a time is not supported.")
            return indices_of_the_stars[0]
        return indices_of_the_stars
    
    def get_density_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying density profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_density_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def set_density_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting density profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        self._check_supplied_values(len(values), number_of_zones)
        self.set_density_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
    
    def get_radius_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying radius profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_radius_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def set_radius_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting radius profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        self._check_supplied_values(len(values), number_of_zones)
        self.set_radius_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
    
    def get_temperature_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying temperature profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_temperature_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def set_temperature_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting temperature profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        self._check_supplied_values(len(values), number_of_zones)
        self.set_temperature_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
    
    def get_mu_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying mean-molecular-weight profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        return self.get_mu_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_names_of_species(self, indices_of_the_stars, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance names")
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        return list(self.get_name_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1) | units.none
        ).value_in(units.string))

    def get_chemical_abundance_profiles(self, indices_of_the_stars, number_of_zones = None, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        grid = numpy.indices((number_of_species, number_of_zones))
        return self.get_mass_fraction_of_species_at_zone(
            [indices_of_the_stars] * number_of_zones * number_of_species, 
            units.none.new_quantity(grid[0].flatten()+1), 
            units.none.new_quantity(grid[1].flatten())
        ).reshape((number_of_species, number_of_zones))
    
    def set_chemical_abundance_profiles(self, indices_of_the_stars, values, number_of_zones = None, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting chemical abundance profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars).number
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars).number
        self._check_supplied_values(len(values), number_of_species, type_string = "chemical species")
        self._check_supplied_values(len(values[0]), number_of_zones)
        grid = numpy.indices((number_of_species, number_of_zones))
        self.set_mass_fraction_of_species_at_zone(
            [indices_of_the_stars] * number_of_zones * number_of_species, 
            units.none.new_quantity(grid[0].flatten()+1), 
            units.none.new_quantity(grid[1].flatten()),
            values.reshape((number_of_species*number_of_zones, ))
        )

