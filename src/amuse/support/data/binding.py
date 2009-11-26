from amuse.support.data import parameters
import numpy

from amuse.support.units import nbody_system

class InterfaceWithParametersBinding(object):
    parameter_definitions = []
    
    def __init__(self, convert_nbody = None):
               
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        

class InterfaceWithObjectsBinding(object):
    def __init__(self):
        self.mapping_from_particleid_to_index = {}
        
    
    def setup_particles(self, particles):
        keyword_arguments = {}
        
        attributes = []
        units = []
        keywords = []
        
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.is_required_for_setup():
                attribute_definition.for_setup_fill_arguments_for_attributelist_get(
                    attributes,
                    units,
                    keywords,
                )
        
        def convert_to_nbody(x):
            if nbody_system.is_nbody_unit(x):
                return self.convert_nbody.unit_to_unit_in_si(x)
            else:
                return x
                
        units = map(convert_to_nbody, units)
        
        list_of_values = particles.attributelist.get_values_of_all_particles_in_units(attributes, units)
    
        
        for keyword, values in zip(keywords, list_of_values):
            keyword_arguments[keyword] = values
            
        indices, errors = self.new_particle(**keyword_arguments)
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not setup a particle")
        
        d = {}
        
        for index in range(len(particles)):
            d[particles[index]] = indices[index]
            
        self.mapping_from_particleid_to_index = d
        

    def update_particles(self, particles):
        ids = map(lambda x : self.mapping_from_particleid_to_index[x], particles)
        
        states = self.get_state(ids)
        
        attributes = []
        units = []
        values = []
        
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.is_required_for_setup():
                attribute_definition.for_state_fill_arguments_for_attributelist_set(
                    states,
                    attributes,
                    units,
                    values,
                )
        
        def convert_to_nbody(x):
            if nbody_system.is_nbody_unit(x):
                return self.convert_nbody.unit_to_unit_in_si(x)
            else:
                return x
          
        
        units = map(convert_to_nbody, units)
        
        time = self.current_model_time()
        
        particles.attributelist.set_values_of_particles_in_units(particles.particles, attributes, values, units, time)
        
    
    def set_attribute(self, attribute_name, particles):
        attribute_definition = self.get_attribute_definition(attribute_name)
        
        if attribute_definition is None:
            return
            
        attributes = []
        units = []
        keywords = []
        
        attribute_definition.for_setter_fill_arguments_for_attributelist_get(
            attributes,
            units,
            keywords,
        )
        
        
        def convert_to_nbody(x):
            if nbody_system.is_nbody_unit(x):
                return self.convert_nbody.unit_to_unit_in_si(x)
            else:
                return x
          
        units = map(convert_to_nbody, units)
        
        list_of_values = particles.attributelist.get_values_of_all_particles_in_units(attributes, units)
        
        keyword_arguments = {}
        for keyword, values in zip(keywords, list_of_values):
            keyword_arguments[keyword] = values
            
        ids = map(lambda x : self.mapping_from_particleid_to_index[x], particles)
        keyword_arguments['id'] = ids
        
        print keyword_arguments
        errors = getattr(self, attribute_definition.setter[0])(**keyword_arguments)

    def update_attribute(self, attribute_name, particles):
        attribute_definition = self.get_attribute_definition(attribute_name)
        
        if attribute_definition is None:
            return
            
        ids = map(lambda x : self.mapping_from_particleid_to_index[x], particles)
        
        results = getattr(self, attribute_definition.getter[0])(ids)
        
        if isinstance(results, tuple):
            results, errors = results
        
        attributes = []
        units = []
        values = []
        
        attribute_definition.for_getter_fill_arguments_for_attributelist_set(
            results,
            attributes,
            units,
            values,
        )
        
        def convert_to_nbody(x):
            if nbody_system.is_nbody_unit(x):
                return self.convert_nbody.unit_to_unit_in_si(x)
            else:
                return x
          
        
        units = map(convert_to_nbody, units)
        
        time = self.current_model_time()
        
        particles.attributelist.set_values_of_particles_in_units(particles.particles, attributes, values, units, time)
                
    
    def get_attribute_definition(self, attribute_name):
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.name == attribute_name:
                return attribute_definition
        return None
            
    def current_model_time(self):
        raise AttributeError("Must implement current_model_time method")
        
