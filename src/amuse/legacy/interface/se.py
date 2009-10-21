"""
Stellar Dynamics Interface Defintion
"""

from amuse.legacy.support.core import legacy_function, LegacyFunctionSpecification

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
    def get_type():
        """
        Retrieve the type of the star. The meaning of the stellar type is defined
        by the code. (Difference between stellar type and type must be explained)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the type of")
        function.addParameter('type', dtype='int32', direction=function.OUT
            , description="The type. ")
        function.result_type = 'i'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A star with the given index was not found.
        """
        return function
        
    @legacy_function   
    def get_stellar_type():
        """
        Retrieve the stellar type of the star. The meaning of the stellar type is defined
        by the code. (Difference between stellar type and type must be explained)
        """
        function = LegacyFunctionSpecification()  
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
    def new_zams_star():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('mass', dtype='float64', direction=function.IN
            , description="The initial mass of the star")
        function.addParameter('age_tag', dtype='float64', direction=function.IN
            , description="Starting age of the star *to be specified exactly*")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            evolution could not complete, solution did not converge
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
    def initialize_stars():
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
    