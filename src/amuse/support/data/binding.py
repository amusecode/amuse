from amuse.support.data import parameters
import numpy

from amuse.support.units import nbody_system

class InterfaceWithParametersBinding(object):
    parameter_definitions = []
    
    def __init__(self, convert_nbody = None):
               
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        


class ParticlesInterface(object):
    
    def __init__(self, code_interface):
        self.code_interface = code_interface
        
    #def __getattr__(self, name_of_the_attribute):
    #    pass
        
    #def __setattr__(self, name_of_the_attribute, value):
    #    pass
        
    def __len__(self):
        return interface.get_len()
    
    def iter_particles(self):
        pass
        
    def add_particles(self, particles):
        attribute_definitions = filter(
            lambda x : x.is_required_for_setup(), 
            self.code_interface.attribute_definitions
        )
        
        attributes = []
        units = []
        keywords = []
        
        for attribute_definition in attribute_definitions:
            attribute_definition.for_setup_fill_arguments_for_attributelist_get(
                attributes,
                units,
                keywords,
            )
        
        units = map(self.code_interface.convert_to_nbody, units)
        
        lists_of_values = particles.attributelist.get_values_of_all_particles_in_units(attributes, units)
        
        keyword_arguments = {}
        for keyword, values in zip(keywords, lists_of_values):
            keyword_arguments[keyword] = values
            
        indices, errors = self.code_interface.new_particle(**keyword_arguments)
        
        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not setup a particle")
        
        
        keys = particles.attributelist.particle_keys
        d = {}
        for index in range(len(particles)):
            d[keys[index]] = indices[index]
            
        self.mapping_from_particle_key_to_index = d
        
    def copy_to(self, particles):
        attribute_definitions = filter(
            lambda x : x.is_required_for_setup(), 
            self.code_interface.attribute_definitions
        )
        
        indices_in_the_code = []
        
        keys = particles.attributelist.particle_keys
        selected_keys = []
        for particle_key in particles.attributelist.particle_keys:
            if particle_key in self.mapping_from_particle_key_to_index:
                indices_in_the_code.append(self.mapping_from_particle_key_to_index[particle_key])
                selected_keys.append(particle_key)
            else:
                pass #raise an exception?
                
        keyword_results = self.code_interface.get_state(indices_in_the_code)
        
        if '__result' in keyword_results:
            errorcodes = keyword_results['__result']
            pass
        
        attributes = []
        units = []
        values = []
        
        for attribute_definition in attribute_definitions:
            attribute_definition.for_state_fill_arguments_for_attributelist_set(
                keyword_results,
                attributes,
                units,
                values,
            )
        
        def convert_to_nbody(x):
            if nbody_system.is_nbody_unit(x):
                return self.convert_nbody.unit_to_unit_in_si(x)
            else:
                return x
          
        
        units = map(self.code_interface.convert_to_nbody, units)
        
        time = self.code_interface.current_model_time()
        
        particles.attributelist.set_values_of_particles_in_units(selected_keys, attributes, values, units, time)
    
    def copy_values_of_attribute_to(self, attribute_name, particles):
        attribute_definition = self.code_interface.get_attribute_definition(attribute_name)
        
        if attribute_definition is None:
            return
            
        indices_in_the_code = []
        
        keys = particles.attributelist.particle_keys
        selected_keys = []
        for particle_key in particles.attributelist.particle_keys:
            if particle_key in self.mapping_from_particle_key_to_index:
                indices_in_the_code.append(self.mapping_from_particle_key_to_index[particle_key])
                selected_keys.append(particle_key)
            else:
                pass #raise an exception?
        
        results = getattr(self.code_interface, attribute_definition.getter[0])(indices_in_the_code)
        
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
        
        units = map(self.code_interface.convert_to_nbody, units)
        
        time = self.code_interface.current_model_time()
        
        particles.attributelist.set_values_of_particles_in_units(selected_keys, attributes, values, units, time)
    
    
    
    def take_values_of_attribute_from(self, attribute_name, particles):
        attribute_definition = self.code_interface.get_attribute_definition(attribute_name)
        
        if attribute_definition is None:
            return
            
        indices_in_the_code = []
        
        keys = particles.attributelist.particle_keys
        selected_keys = []
        for particle_key in particles.attributelist.particle_keys:
            if particle_key in self.mapping_from_particle_key_to_index:
                indices_in_the_code.append(self.mapping_from_particle_key_to_index[particle_key])
                selected_keys.append(particle_key)
            else:
                pass #raise an exception?
                
        attributes = []
        units = []
        keywords = []
        
        attribute_definition.for_setter_fill_arguments_for_attributelist_get(
            attributes,
            units,
            keywords,
        )
        
        units = map(self.code_interface.convert_to_nbody, units)
        
        list_of_values = particles.attributelist.get_values_of_all_particles_in_units(attributes, units)
        
        keyword_arguments = {}
        for keyword, values in zip(keywords, list_of_values):
            keyword_arguments[keyword] = values
            
        
        keyword_arguments['id'] = indices_in_the_code
        keyword_arguments['index_of_the_particle'] = indices_in_the_code
        
        setter_method = getattr(self.code_interface, attribute_definition.setter[0])
        
        if len(indices_in_the_code) > 0: 
            if hasattr(setter_method, "specification"):
                if not setter_method.specification.can_handle_array:
                    raise Exception(
                        "setter method <{0}> of a '{1}' object, cannot handle arrays".format(attribute_definition.setter[0], type(self).__name__)
                    )
                    
        errors = setter_method(**keyword_arguments)
                
        
        
class InterfaceWithObjectsBinding(object):
    def __init__(self):
        self.mapping_from_particleid_to_index = {}
        self.particles = ParticlesInterface(self)
    
    def convert_to_nbody(self, x):
        if nbody_system.is_nbody_unit(x):
            return self.convert_nbody.unit_to_unit_in_si(x)
        else:
            return x
                
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def update_particles(self, particles):
        self.particles.copy_to(particles)
    
    def set_attribute(self, attribute_name, particles):
        self.particles.take_values_of_attribute_from(attribute_name, particles)

    def update_attribute(self, attribute_name, particles):
        self.particles.copy_values_of_attribute_to(attribute_name, particles) 
    
    def get_attribute_definition(self, attribute_name):
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.name == attribute_name:
                return attribute_definition
        return None
            
    def current_model_time(self):
        raise AttributeError("Must implement current_model_time method")
        
