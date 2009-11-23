from amuse.support.data import parameters
import numpy

class InterfaceWithParametersBinding(object):
    parameter_definitions = []
    
    def __init__(self, convert_nbody = None):
               
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        

class InterfaceWithObjectsBinding(object):
   
    
    def setup_particles(self, particles):
        keyword_arguments = {}
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.is_required_for_setup():
                values = particles.get_values_of_attribute(attribute_definition.name)
                attribute_definition.set_keyword_arguments(self, values, keyword_arguments)
        
        print keyword_arguments
        indices, errors = self.new_particle(**keyword_arguments)
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not setup a particle")
        
        module_id = id(self)
        for particle, new_index in zip(particles, indices):
            particle._module_ids_to_index[module_id] = (self, new_index)
            

    def update_particles(self, particles):
        ids = list(particles.ids_for_module_with_id(id(self)))
        state = self.get_state(ids)
        time = self.current_model_time()
        for attribute_definition in self.attribute_definitions:
            values = attribute_definition.get_keyword_results(self, state)
            particles.set_values_of_attribute(attribute_definition.name, time, values)
            
    
    def set_attribute(self, attribute_name, particles):
        attribute_definition = self.get_attribute_definition(attribute_name)
        
        if attribute_definition is None:
            return
            
        values = particles.get_values_of_attribute_for_module_with_id(
            id(self), 
            attribute_definition.name
        )
        ids = list(particles.ids_for_module_with_id(id(self)))
        
        attribute_definition.set_values(self, ids, values)

    def update_attribute(self, attribute_name, particles):
        attribute_definition = self.get_attribute_definition(attribute_name)
        
        if attribute_definition is None:
            return
            
        time = self.current_model_time()
        
        ids =  list(particles.ids_for_module_with_id(id(self)))
        values = attribute_definition.get_values(self, ids)
        
        particles.set_values_of_attribute_for_module_with_id(
            id(self), 
            attribute_definition.name, 
            time, 
            values)
                
    
    def get_attribute_definition(self, attribute_name):
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.name == attribute_name:
                return attribute_definition
        return None
            
    def current_model_time(self):
        raise AttributeError("Must implement current_model_time method")
        
