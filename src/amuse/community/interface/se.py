"""
Stellar Dynamics Interface Defintion
"""

import numpy

from amuse.units import units
from amuse.datamodel import Particles, Particle
from amuse.support import exceptions
from amuse.ext.spherical_model import EnclosedMassInterpolator
from amuse.community.interface import common

from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification
class StellarEvolutionInterface(common.CommonCodeInterface):


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
    def evolve_one_step():
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
    def evolve_for():
        """
        Evolve the star for exactly the given time period.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to evolve")
        function.addParameter('delta_t', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
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
    
    @legacy_function
    def recommit_particles():
        """
        Let the code perform reinitialization actions after additional particles 
        have been added or removed. Not every code needs this functionality, and 
        it may be possible to create stars after this function has been called.
        """
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function  
    


class StellarEvolution(common.CommonCode):
    
    def evolve_model(self, end_time = None, keep_synchronous = True):
        if not keep_synchronous:
            for particle in self.particles:
                particle.evolve_one_step()
            return
        
        delta_time = end_time-self.model_time if end_time else 0.99*min(self.particles.time_step)
        for particle in self.particles:
            particle.evolve_for(delta_time)
        self.model_time += delta_time
    
    def _evolve_model_old(self, end_time = None, keep_synchronous = True):
        """
        This is the old implementation of evolve_model. Even with (keep_synchronous = True) 
        it is unable to evolve all stars to a common age, since it relies on the 
        individual timesteps as determined by the community code. Furthermore, it 
        is not suited to simulations with ongoing star formation, since it evolves 
        newly created stars to the same age as the old stars. Finally, this old 
        implementation has substantially more communication overhead. 
        """
        if end_time is None:
            if keep_synchronous:
                ages = self.particles.age
                index, min_age = min(enumerate(ages), key=itemgetter(1))
                self.particles[index].evolve_one_step()
                new_age = self.particles[index].age
                for particle in self.particles.select(lambda x : x < new_age, ["age"]):
                    particle.evolve_one_step()
            else:
                self.particles.evolve_one_step()
        else:
            for particle in self.particles:
                while particle.age < end_time:
                    particle.evolve_one_step()
    
    def define_state(self, object):
        common.CommonCode.define_state(self, object)
        
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        object.add_method('UPDATE','before_get_parameter')
        object.add_method('EVOLVED','before_get_parameter')
        
        
        
        object.add_method('EDIT', 'new_particle')
        object.add_method('EDIT', 'delete_star')
        object.add_method('UPDATE', 'new_particle')
        object.add_method('UPDATE', 'delete_star')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_particle', False)
        object.add_transition('RUN', 'UPDATE', 'finalize_stellar_model', False)
        object.add_transition('RUN', 'UPDATE', 'delete_star', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        object.add_method('RUN', 'evolve_model')
        object.add_method('RUN', 'evolve_for')
        object.add_method('RUN', 'evolve_one_step')
        object.add_method('RUN', 'get_age')
        object.add_method('RUN', 'get_mass')
        object.add_method('RUN', 'get_luminosity')
        object.add_method('RUN', 'get_radius')
        object.add_method('RUN', 'get_stellar_type')
        object.add_method('RUN', 'get_temperature')
    
    def define_methods(self, object):
        common.CommonCode.define_methods(self, object)
        object.add_method(
            "evolve_one_step",
            (object.INDEX,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "evolve_for",
            (object.INDEX, units.yr),
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
            "get_metallicity", 
            (), 
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_metallicity", 
            (object.NO_UNIT, ), 
            (object.ERROR_CODE,)
        )


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
        object.add_method(set_name, 'calculate_core_mass')
        object.add_method(set_name, 'calculate_helium_exhausted_core_mass')
        object.add_getter(set_name, 'get_central_temperature', names = ('central_temperature',))
        object.add_getter(set_name, 'get_central_density', names = ('central_density',))
    
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
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "get_temperature_at_zone", 
            (object.INDEX,object.NO_UNIT,), 
            (units.K, object.ERROR_CODE,)
        )
        object.add_method(
            "set_temperature_at_zone", 
            (object.INDEX, object.NO_UNIT, units.K,), 
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_density_at_zone", 
            (object.INDEX,object.NO_UNIT,), 
            (units.g/units.cm**3, object.ERROR_CODE,)
        )
        object.add_method(
            "set_density_at_zone", 
            (object.INDEX, object.NO_UNIT, units.g/units.cm**3,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_radius_at_zone", 
            (object.INDEX,object.NO_UNIT,), 
            (units.cm, object.ERROR_CODE,)
        )
        object.add_method(
            "set_radius_at_zone", 
            (object.INDEX, object.NO_UNIT, units.cm,), 
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_mu_at_zone", 
            (object.INDEX, object.NO_UNIT,), 
            (units.amu, object.ERROR_CODE,)
        )
        object.add_method(
            "get_number_of_species", 
            (object.INDEX,), 
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "get_name_of_species", 
            (object.INDEX, object.NO_UNIT,), 
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "get_mass_fraction_of_species_at_zone", 
            (object.INDEX, object.NO_UNIT, object.NO_UNIT,), 
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_mass_fraction_of_species_at_zone", 
            (object.INDEX, object.NO_UNIT, object.NO_UNIT, object.NO_UNIT,), 
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
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_density_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_central_density(self, indices_of_the_stars):
        return self.get_density_at_zone(indices_of_the_stars, 0)
    
    def set_density_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting density profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_density_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
        if hasattr(self, "_erase_memory"):
            self._erase_memory(indices_of_the_stars)
    
    def get_radius_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying radius profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_radius_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def set_radius_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting radius profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_radius_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
        if hasattr(self, "_erase_memory"):
            self._erase_memory(indices_of_the_stars)
    
    def get_temperature_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying temperature profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_temperature_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_central_temperature(self, indices_of_the_stars):
        return self.get_temperature_at_zone(indices_of_the_stars, 0)
    
    def set_temperature_profile(self, indices_of_the_stars, values, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting temperature profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_zones)
        self.set_temperature_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none, values)
        if hasattr(self, "_erase_memory"):
            self._erase_memory(indices_of_the_stars)
    
    def get_mu_profile(self, indices_of_the_stars, number_of_zones = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying mean-molecular-weight profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        return self.get_mu_at_zone([indices_of_the_stars]*number_of_zones, range(number_of_zones) | units.none)
    
    def get_names_of_species(self, indices_of_the_stars, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance names")
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        return list(self.get_name_of_species(
            [indices_of_the_stars]*number_of_species, 
            range(1,number_of_species+1) | units.none
        ))

    def get_chemical_abundance_profiles(self, indices_of_the_stars, number_of_zones = None, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying chemical abundance profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        grid = numpy.indices((number_of_species, number_of_zones))
        return self.get_mass_fraction_of_species_at_zone(
            [indices_of_the_stars] * number_of_zones * number_of_species, 
            grid[0].flatten()+1, 
            grid[1].flatten()
        ).reshape((number_of_species, number_of_zones))
    
    def set_chemical_abundance_profiles(self, indices_of_the_stars, values, number_of_zones = None, number_of_species = None):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Setting chemical abundance profiles")
        if number_of_zones is None:
            number_of_zones = self.get_number_of_zones(indices_of_the_stars)
        if number_of_species is None:
            number_of_species = self.get_number_of_species(indices_of_the_stars)
        self._check_supplied_values(len(values), number_of_species, type_string = "chemical species")
        self._check_supplied_values(len(values[0]), number_of_zones)
        grid = numpy.indices((number_of_species, number_of_zones))
        self.set_mass_fraction_of_species_at_zone(
            [indices_of_the_stars] * number_of_zones * number_of_species, 
            grid[0].flatten()+1, 
            grid[1].flatten(),
            values.reshape((number_of_species*number_of_zones, ))
        )
        if hasattr(self, "_erase_memory"):
            self._erase_memory(indices_of_the_stars)
    
    def calculate_core_mass(self, indices_of_the_stars, species=None, core_H_abundance_limit=1.0e-4, split_species=False):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying the core mass")
        chemical_abundance_profiles = self.get_chemical_abundance_profiles(indices_of_the_stars)
        index_core = numpy.searchsorted(chemical_abundance_profiles[0], core_H_abundance_limit)
        densities, radii_cubed = self._get_densities_radii_cubed(indices_of_the_stars, index_core, 
            split_species, species, chemical_abundance_profiles)
        sum_axis = 1 if split_species else None
        return (numpy.pi * 4.0/3.0 * (densities * (radii_cubed[1:] - radii_cubed[:-1])).sum(axis=sum_axis)).as_quantity_in(units.MSun)
    
    def _get_index_helium_exhausted_core(self, indices_of_the_stars, chemical_abundance_profiles, core_He_abundance_limit):
        helium_abundance_profile = 0 * chemical_abundance_profiles[0]
        for i, species_name in enumerate(self.get_names_of_species(indices_of_the_stars)):
            if "he" in species_name or "He" in species_name:
                helium_abundance_profile += chemical_abundance_profiles[i]
        return numpy.searchsorted(helium_abundance_profile, core_He_abundance_limit)
    
    def _get_densities_radii_cubed(self, indices_of_the_stars, index_core, 
            split_species, species, chemical_abundance_profiles):
        densities = self.get_density_profile(indices_of_the_stars)[:index_core]
        if split_species:
            if species is None:
                species_profiles = chemical_abundance_profiles[..., :index_core]
            else:
                species_profiles = []
                for i, species_name in enumerate(self.get_names_of_species(indices_of_the_stars)):
                    if species_name in species:
                        species_profiles.append(chemical_abundance_profiles[i, :index_core])
            densities = densities * species_profiles
        else:
            if not species is None:
                fraction = 0 * chemical_abundance_profiles[0, :index_core]
                for i, species_name in enumerate(self.get_names_of_species(indices_of_the_stars)):
                    if species_name in species:
                        fraction += chemical_abundance_profiles[i, :index_core]
                densities = densities * fraction
        
        radii = self.get_radius_profile(indices_of_the_stars)[:index_core]
        radii.prepend(0 | units.m)  
        radii_cubed = radii**3
        return densities, radii_cubed
    
    def calculate_helium_exhausted_core_mass(self, indices_of_the_stars, species=None, core_He_abundance_limit=1.0e-4, split_species=False):
        indices_of_the_stars = self._check_number_of_indices(indices_of_the_stars, action_string = "Querying the core mass")
        chemical_abundance_profiles = self.get_chemical_abundance_profiles(indices_of_the_stars)
        index_core = self._get_index_helium_exhausted_core(indices_of_the_stars, chemical_abundance_profiles, core_He_abundance_limit)
        densities, radii_cubed = self._get_densities_radii_cubed(indices_of_the_stars, index_core, 
            split_species, species, chemical_abundance_profiles)
        sum_axis = 1 if split_species else None
        return (numpy.pi * 4.0/3.0 * (densities * (radii_cubed[1:] - radii_cubed[:-1])).sum(axis=sum_axis)).as_quantity_in(units.MSun)
    
    def merge_colliding(self, primaries, secondaries, collision_code, 
            code_options=dict(), code_parameters=dict(), return_merge_products=["se", "gd"], create_new_key=True):
        return merge_colliding_in_stellar_evolution_code(self, 
            primaries, secondaries, collision_code, 
            code_options=code_options, code_parameters=code_parameters, 
            return_merge_products=return_merge_products, create_new_key=create_new_key)
    

def merge_colliding_in_stellar_evolution_code(stellar_evolution_code, primaries, secondaries, collision_code, 
        code_options=dict(), code_parameters=dict(), return_merge_products=["se", "gd"], create_new_key=True):
    primaries = primaries.as_set()
    secondaries = secondaries.as_set()
    star_collider = collision_code(**code_options)
    for (par_name, value) in code_parameters.iteritems():
        setattr(star_collider.parameters, par_name, value)
    star_collider.commit_parameters()
    star_collider.particles.add_particles(primaries)
    star_collider.particles.add_particles(secondaries)
    se_colliders = star_collider.particles.get_intersecting_subset_in(stellar_evolution_code.particles)
    for col_particle, se_particle in zip(star_collider.particles, se_colliders):
        number_of_zones     = se_particle.get_number_of_zones()
        mm1                 = se_particle.get_mass_profile(number_of_zones = number_of_zones)* se_particle.mass
        mass_profile        = se_particle.get_cumulative_mass_profile(number_of_zones = number_of_zones) * se_particle.mass
        density_profile     = se_particle.get_density_profile(number_of_zones = number_of_zones)
        radius_profile      = se_particle.get_radius_profile(number_of_zones = number_of_zones)
        temperature_profile = se_particle.get_temperature_profile(number_of_zones = number_of_zones)
        lum                 = se_particle.get_luminosity_profile(number_of_zones = number_of_zones)
        pressure_profile    = se_particle.get_pressure_profile(number_of_zones = number_of_zones)
        mu_profile          = se_particle.get_mu_profile(number_of_zones = number_of_zones)
        composition_profile = se_particle.get_chemical_abundance_profiles(number_of_zones = number_of_zones)
        col_particle.add_shell(mm1,mass_profile, radius_profile, density_profile, 
            pressure_profile, temperature_profile,lum, mu_profile, composition_profile[0], 
            composition_profile[1]+composition_profile[2], composition_profile[3], 
            composition_profile[4], composition_profile[5], composition_profile[6], 
            composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
    
    stellar_evolution_code.particles.remove_particles(star_collider.native_stars)
    
    gd_merge_products = Particles()
    for primary, secondary in zip(primaries, secondaries):
        merge_product = Particle()
        merge_product.primary = primary
        merge_product.secondary = secondary
        new_particle = star_collider.merge_products.add_particle(merge_product)
        stellar_model = new_particle.internal_structure()
        if create_new_key:
            new_key = merge_product.key
        else:
            new_key = primary.key if primary.mass > secondary.mass else secondary.key
        stellar_evolution_code.new_particle_from_model(stellar_model, 0.0|units.Myr, key=new_key)
        
        if "gd" in return_merge_products:
            merge_product = Particle(key=new_key)
            merge_product.mass = stellar_model.mass[-1] 
            merge_product.radius= stellar_model.radius[-1]
            gd_colliders = (primary + secondary)
            merge_product.position = gd_colliders.center_of_mass()
            merge_product.velocity = gd_colliders.center_of_mass_velocity()
            gd_merge_products.add_particle(merge_product)
        
    star_collider.stop()
    
    result = []
    for type in return_merge_products:
        if type == "se":
            result.append(gd_merge_products.get_intersecting_subset_in(stellar_evolution_code.particles))
        elif type == "gd":
            result.append(gd_merge_products)
        else:
            print "Unexpected value in return_merge_products, must be 'gd' (gravity particles) or 'se' (stellar evolution particles):", type
    return result

