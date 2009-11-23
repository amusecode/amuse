from amuse.legacy import *
from amuse.legacy.interface.se import StellarEvolution
from amuse.legacy.support.lit import LiteratureRefs
from amuse.support.data.binding import InterfaceWithParametersBinding, InterfaceWithObjectsBinding

import os

class EVtwin(LegacyInterface, LiteratureRefs, StellarEvolution):
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
    def get_maximum_number_of_stars():
        """
        Retrieve the maximum number of stars that can be
        handled by this instance.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('metallicity', dtype='int32', direction=function.OUT,
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
        function.addParameter('metallicity', dtype='int32', direction=function.IN,
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
        
class EVtwinBinding(InterfaceWithParametersBinding, InterfaceWithObjectsBinding):
    
    def __init__(self):
        InterfaceWithParametersBinding.__init__(self)
        InterfaceWithObjectsBinding.__init__(self)
        
    parameter_definitions = [
        parameters.ModuleMethodParameterDefinition_Next(
            "get_maximum_number_of_stars",
            "set_maximum_number_of_stars",
            "maximum_number_of_stars", 
            "Maximum number of stars that can be allocated", 
            units.none, 
            10 | units.none
        ),
        
        parameters.ModuleMethodParameterDefinition_Next(
            "get_metallicity",
            "set_metallicity",
            "metallicity", 
            "Metallicity of all stats", 
            units.percentage, 
            0.02 | units.percentage
        ),
        
        
        parameters.ModuleMethodParameterDefinition_Next(
            None,
            "set_ev_path",
            "path_to_data", 
            "Path to the data directory", 
            units.string, 
            "src" | units.string
        ),
        
    ]
    
    attribute_definitions = [
        attributes.ScalarAttributeDefinition_Next(
            None,
            "get_stellar_type",
            None,
            "type",
            "star type",
             units.stellar_type,
             1 | units.stellar_type
        ),
        attributes.ScalarAttributeDefinition_Next(
            None,
            "get_mass",
            "mass",
            "mass",
            "mass of the star",
             units.MSun,
             0.0 | units.MSun
        ),
        attributes.ScalarAttributeDefinition_Next(
            None,
            "get_age",
            None,
            "age",
            "current age of the star",
             units.Myr,
             1.0 | units.Myr
        ),
        attributes.ScalarAttributeDefinition_Next(
            None,
            "get_radius",
            None,
            "radius",
            "current radius of the star",
             units.RSun,
             1.0 | units.RSun
        ),
        attributes.ScalarAttributeDefinition_Next(
            None,
            "get_luminosity",
            None,
            "luminosity",
            "current luminosity of the star",
             units.LSun,
             1.0 | units.LSun
        ),
    ]
    
    
    def update_particles(self, particles):
        ids = list(particles.ids_for_module_with_id(id(self)))
        attribute_definition = self.get_attribute_definition("age");
        times =  attribute_definition.get_values(self, ids)
        
        module_id = id(self)
        for attribute_definition in self.attribute_definitions:
            values = attribute_definition.get_values(self, ids)
            particles.set_values_of_attribute_for_module_with_id(module_id, attribute_definition.name, values=values, times=times)
            
    def evolve_particles(self, particles, end_time):
        module_id = id(self)
        
        
        for particle in particles:
            if not module_id in particle._module_ids_to_index:
                    continue
                    
            index = particle._module_ids_to_index[module_id][1]
            
                
            current_age, error = self.get_age(index)
            current_age |= units.Myr
            
            print index, current_age
            
            while current_age < end_time:
                
                errorcode = self.evolve(index)
                if errorcode < 0:
                    raise Exception("Error during evolution of a star, errorcode = " + errorcode)
                
                current_age, error = self.get_age(index)
                current_age |= units.Myr
                
                print current_age
                
                for attribute_defintion in self.attribute_definitions:
                    values = attribute_definition.get_values(self, [index])
                    particle.set_value_of_attribute(
                        attribute_definition.name,
                        values[0]
                    )
                
        
