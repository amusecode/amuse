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
        attributes.AttributeDefinition(
            name = "type",
            getter = ("get_stellar_type", ["stellar_type"]),
            description = "star type",
            unit = units.stellar_type,
            default = 1 | units.stellar_type             
        ),
        attributes.AttributeDefinition(
            name = "mass",
            setup_parameters = ["mass"],
            getter = ("get_mass", ["mass"]),
            description = "mass of the star",
            unit = units.MSun,
            default = 0.0 | units.MSun             
        ),
        attributes.AttributeDefinition(
            name = "age",
            getter = ("get_age", ["age"]),
            description = "current age of the star",
            unit = units.Myr,
            default = 1.0 | units.Myr ,            
        ),
        attributes.AttributeDefinition(
            name = "radius",
            getter = ("get_radius", ["radius"]),
            description = "current radius of the star",
            unit = units.RSun,
            default = 1.0 | units.RSun             
        ),
        attributes.AttributeDefinition(
            name = "luminosity",
            getter = ("get_luminosity", ["luminosity"]),
            description = "current luminosity of the star",
            unit = units.LSun,
            default = 1.0 | units.LSun             
        ),
    ]
    
    
    def update_particles(self, particles):
        self._current_model_time = None
        
        for attribute_definition in self.attribute_definitions:
            self.update_attribute(attribute_definition.name, particles)
            
    def evolve_particles(self, particles, end_time):
        module_id = id(self)
        
        
        for particle in particles:
            if not module_id in particle._module_ids_to_index:
                    continue
                    
            index = particle._module_ids_to_index[module_id][1]
            
                
            current_age, error = self.get_age(index)
            current_age |= units.Myr
            
            while current_age < end_time:
                
                errorcode = self.evolve(index)
                if errorcode < 0:
                    raise Exception("Error during evolution of a star, errorcode = " + errorcode)
                
                current_age, error = self.get_age(index)
                current_age |= units.Myr
                
                for attribute_defintion in self.attribute_definitions:
                    values = attribute_definition.get_values(self, [index])
                    particle.set_value_of_attribute(
                        attribute_definition.name,
                        values[0]
                    )
                    
    def current_model_time(self):
        return self._current_model_time
                
    def get_state(self, indices):
        return self.get_mass(indices)
