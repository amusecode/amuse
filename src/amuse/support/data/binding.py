from amuse.support.data import parameters
from amuse.support.data.core import Particles, ParticleInformationChannel, Particle
import numpy

from amuse.support.units import nbody_system

class InterfaceWithParametersBinding(object):
    parameter_definitions = []
    
    def __init__(self, convert_nbody = None):
               
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        

class ParticlesInTheModule(object):
    
    def __init__(self, code_interface):
        self.code_interface = code_interface
        self.mapping_from_particle_key_to_index_in_the_code = {}
        self.particle_keys = []
                
    def __len__(self):
        return interface.get_len()
        
    def complies_with_state(self, attributes):
        attributes_in_getstate = self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
        print attributes_in_getstate
        if not attributes_in_getstate:
            return False
        
        if len(attributes_in_getstate) == len(attributes):
            check = set(attributes)
            for x in attributes_in_getstate:
                if x in check:
                    check.remove(x)
                else:
                    return False
            return len(check) == 0
        else:
            return False
            
    def select_getters_for(self, attributes):
        attributes_in_getstate = self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
        if self.complies_with_state(attributes):
            return [("get_state", self._mapping_from_attribute_names_to_set_keyword_and_unit())]
        set_of_attributes = set(attributes)
        
        result = []
        for x in self.code_interface.attribute_definitions:
            y = x.attributes_in_getter(set_of_attributes)
            if not y is None:
                result.append(y)
                
        if set_of_attributes:
            raise Exception("Do not have attributes {0}".format(sorted(set_of_attributes)))
        return result
        
    
    def select_setters_for(self, attributes):
        attributes_in_getstate = self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
        if self.complies_with_state(attributes):
            return [("set_state", self._mapping_from_attribute_names_to_set_keyword_and_unit())]
        set_of_attributes = set(attributes)
        
        result = []
        for x in self.code_interface.attribute_definitions:
            y = x.attributes_in_setter(set_of_attributes)
            if not y is None:
                result.append(y)
                
        if set_of_attributes:
            raise Exception("Cannot set attributes {0}".format(sorted(set_of_attributes)))
        return result
        
    def _mapping_from_attribute_names_to_set_keyword_and_unit(self):
        result = {}
        for attribute_definition in self.code_interface.attribute_definitions:
            if not attribute_definition.is_required_for_setup():
                continue
            result.update(attribute_definition.state_mapping_from_name_to_keyword_and_unit())
        return result
            
    def _set_particles(self, keys, attributes = [], values = []):
        self.particle_keys = keys
        
        mapping = self._mapping_from_attribute_names_to_set_keyword_and_unit()
        keyword_arguments = {}
        for attribute, quantity in zip(attributes, values):
            if attribute in mapping:
                keyword, unit = mapping[attribute]
                unit = self.code_interface.convert_to_nbody(unit)
                keyword_arguments[keyword] = quantity.value_in(unit)
            
        indices, errors = self.code_interface.new_particle(**keyword_arguments)

        for errorcode in errors:
            if errorcode < 0:
                raise Exception("Could not setup a particle")

        d = self.mapping_from_particle_key_to_index_in_the_code
        index = 0
        for key in keys:
            d[key] = indices[index]
            index = index + 1
        
    def _get_values(self, keys, attributes):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
        
        mapping_from_attribute_to_result = {}
        
        for method_name, mapping in self.select_getters_for(attributes):
            keyword_to_attribute_and_unit_mapping = {}
            
            
            method = getattr(self.code_interface, method_name)
            
            if len(indices_in_the_code) > 1: 
                if hasattr(method, "specification") and not method.specification.can_handle_array:
                    raise Exception(
                        "getter method <{0}> of a '{1}' object, cannot handle arrays".format(method_name, type(self).__name__)
                    )
            
            keyword_results = method(indices_in_the_code)
            
            for attribute in mapping.keys():
                keyword, unit = mapping[attribute]
                unit = self.code_interface.convert_to_nbody(unit)
                mapping_from_attribute_to_result[attribute] = unit.new_quantity(keyword_results[keyword])
        
        results = []
        for attribute in attributes:
            results.append(mapping_from_attribute_to_result[attribute])
        return results
        
                
        
    def _set_values(self, keys, attributes, values):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
            
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
            
        for method_name, mapping in self.select_setters_for(attributes):
            keyword_arguments = {}
            
            for attribute, quantity in zip(attributes, values):
                keyword, unit = mapping[attribute]
                unit = self.code_interface.convert_to_nbody(unit)
                keyword_arguments[keyword] = quantity.value_in(unit)
        
            
            method = getattr(self.code_interface, method_name)
            
            if len(indices_in_the_code) > 1: 
                if hasattr(method, "specification") and not method.specification.can_handle_array:
                    raise Exception(
                        "setter method <{0}> of a '{1}' object, cannot handle arrays".format(method_name, type(self.code_interface).__name__)
                    )
                    
            keyword_arguments['id'] = indices_in_the_code
            keyword_arguments['index_of_the_particle'] = indices_in_the_code
            
            errors = method(**keyword_arguments)
        
    def _get_attributes(self):
        return self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()
    
    def _get_keys(self):
        return self.particle_keys

    def _has_key(self, key):
        return key in self.mapping_from_particle_key_to_index_in_the_code
    
    def _state_attributes(self):
        return self._mapping_from_attribute_names_to_set_keyword_and_unit().keys()

    
class InterfaceWithObjectsBinding(object):
    def __init__(self):
        self.mapping_from_particleid_to_index = {}
        self.particles = Particles()
        self.particles.attributelist = ParticlesInTheModule(self)
    
    def convert_to_nbody(self, x):
        if nbody_system.is_nbody_unit(x):
            return self.convert_nbody.unit_to_unit_in_si(x)
        else:
            return x
                
    def setup_particles(self, particles):
        self.particles.add_particles(particles)
        
    def update_particles(self, particles):
        self.particles.copy_values_of_state_attributes_to(particles)
    
    def set_attribute(self, attribute_name, particles):
        particles.copy_values_of_attribute_to(attribute_name, self.particles)
        
    def update_attribute(self, attribute_name, particles):
        self.particles.copy_values_of_attribute_to(attribute_name, particles) 
    
    def get_attribute_definition(self, attribute_name):
        for attribute_definition in self.attribute_definitions:
            if attribute_definition.name == attribute_name:
                return attribute_definition
        return None
            
    def current_model_time(self):
        raise AttributeError("Must implement current_model_time method")
        
