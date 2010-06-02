from amuse.support.data import values
from amuse.support.data.values import Quantity, new_quantity, zero
from amuse.support.units import constants
from amuse.support.units import units
from amuse.support.core import CompositeDictionary


from amuse.support.data.base import *
        
class InMemoryAttributeStorage(AttributeStorage):
    
    def __init__(self):
        self.model_times = None
        self.mapping_from_attribute_to_values_and_unit = {}
        self.mapping_from_particle_to_index = {}
        self.particle_keys = []

    def _set_particles(self, keys, attributes = [], values = []):
        if len(values) != len(attributes):
            raise Exception(
                "you need to provide the same number of value list as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(values)
                )
            )
        if len(values) > 0 and len(keys) != len(values[0]):
            raise Exception(
                "you need to provide the same number of values as particles, found {0} values and {1} particles".format(
                    len(values[0]), len(keys)
                )
            )
        
        if len(self.particle_keys) > 0:
            self.append_to_storage(keys, attributes, values)
        else:
            self.setup_storage(keys, attributes, values)
            
    def setup_storage(self, keys, attributes, values):
        self.mapping_from_attribute_to_values_and_unit = {}
        for attribute, quantity in zip(attributes, values):
            self.mapping_from_attribute_to_values_and_unit[attribute] = AttributeValues(
                attribute,
                quantity.unit,
                quantity.number,
                None
            )
         
        self.particle_keys = numpy.array(keys, dtype='uint64')

        self.reindex()
    
    def append_to_storage(self, keys, attributes, values):
        for attribute, values_to_set in zip(attributes, values):
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                    attribute,
                    values_to_set.unit,
                    length = len(self.particle_keys)
                )
            
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
            values_in_the_right_units = values_to_set.value_in(attribute_values.unit)
            attribute_values.values = numpy.concatenate((attribute_values.values, values_in_the_right_units))
        
        old_length = len(self.particle_keys)
        for attribute_values in self.mapping_from_attribute_to_values_and_unit.values():
            zeros_for_concatenation = numpy.zeros(len(keys), dtype=attribute_values.unit.dtype)
            if len(attribute_values.values) == old_length:
                attribute_values.values = numpy.concatenate((attribute_values.values, zeros_for_concatenation))
        
                
        index = len(self.particle_keys)

        self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys), dtype='uint64')))

        for particle_key in keys:
            self.mapping_from_particle_to_index[particle_key] = index
            index += 1
            
    def _get_values(self, particles, attributes):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute in attributes:
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             if indices is None:
                 selected_values = attribute_values.values
             else:
                 selected_values = attribute_values.values.take(indices)
             results.append(attribute_values.unit.new_quantity(selected_values))
        
        return results
        
    def _set_values(self, particles, attributes, list_of_values_to_set, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, len(indices))
        
        previous_model_times = None
        if list_of_values_to_set is None:
            for attribute in attributes:
                if attribute in self.mapping_from_attribute_to_values_and_unit:
                    attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
                else:
                    raise Exception("unknown attribute '{0}'".format(attribute))
                     
                selected_values = numpy.zeros(len(indices), dtype=attribute_values.values.dtype)
                
                attribute_values.values.put(indices, selected_values)
            return
            
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                   attribute,
                   values_to_set.unit,
                   length = len(self.particle_keys)
                )
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                 
            selected_values = values_to_set.value_in(attribute_values.unit)
            attribute_values.values.put(indices, selected_values)
            if not model_times is None:
                if not previous_model_times is attribute_values.model_times:
                    attribute_values.model_times.put(indices, model_times)
                    previous_model_times = attribute_values.model_times
            
    
    
    def _get_attributes(self):
        return sorted(self.mapping_from_attribute_to_values_and_unit.keys())
    
    
    def _has_key(self, key):
        return key in self.mapping_from_particle_to_index
        
    def _get_keys(self):
        return self.particle_keys
        
    def __len__(self):
        return len(self.particle_keys)
        
    def copy(self):
        copy = InMemoryAttributeStorage()
        copy.mapping_from_particle_to_index = self.mapping_from_particle_to_index.copy()
        copy.particle_keys = self.particle_keys.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_values_and_unit.iteritems():
            copy.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle_key, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            raise AttributeError("particle does not have a "+attribute)
        
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        index = self.mapping_from_particle_to_index[particle_key]
        
        return attribute_values.unit.new_quantity(attribute_values.values[index])
        
            
    def iter_values_of_particle(self, particle_key):
        index = self.mapping_from_particle_to_index[particle_key]
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            yield attribute, (attribute_values.values[index] | attribute_values.unit)
    
    
            
    def iter_values_of(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            return
            
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        values = attribute_values.values
        unit = attribute_values.unit
        particles = self.particle_keys
        
        for index in range(len(self.particle_keys)):
            yield particles[i], (values[i] | unit)
            
   
    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
            
        mapping_from_particle_to_index = self.mapping_from_particle_to_index 
        result = numpy.zeros(len(particles),dtype='int32')
        #result = [mapping_from_particle_to_index[id(particle)] for particle in particles]
        
        index = 0
        for index, particle_key in enumerate(particles):
            result[index] = mapping_from_particle_to_index[particle_key]
            index += 1
        return result
    
    def get_values_of_particles_in_units(self, particles, attributes, target_units):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute, target_unit in zip(attributes, target_units):
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_unit_in_target_unit = attribute_values.unit.value_in(target_unit )
             if indices is None:
                 selected_values = attribute_values.values
             else:
                 selected_values = attribute_values.values.take(indices)
             if value_of_unit_in_target_unit != 1.0:
                selected_values *= value_of_unit_in_target_unit
             results.append(selected_values)
        
        return results
        
    
            
    def _convert_model_times(self, model_times, length):
        if not model_times is None and isinstance(model_times, values.ScalarQuantity):
            return model_times.unit.new_quantity(numpy.linspace(model_times.number, model_times.number, length) )
        else:
            return model_times
    
    def set_values_of_particles_in_units(self, particles, attributes, list_of_values_to_set, source_units, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, len(indices))
        
        previous_model_times = None
        results = []
        for attribute, values_to_set, source_unit in zip(attributes, list_of_values_to_set, source_units):
            selected_values = values_to_set
            
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                   attribute,
                   source_unit,
                   length = len(self.particle_keys)
                )
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                 
            value_of_source_unit_in_list_unit = source_unit.value_in(attribute_values.unit)
            if value_of_source_unit_in_list_unit != 1.0:
                selected_values *= value_of_source_unit_in_list_unit 
             
            attribute_values.values.put(indices, selected_values)
            if not model_times is None:
                if not previous_model_times is attribute_values.model_times:
                    attribute_values.model_times.put(indices, model_times)
                    previous_model_times = attribute_values.model_times
            
        return results
        
             
    def merge_into(self, others):
        source_attributes = []
        source_units = []
        source_valeus = []
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            source_attributes.append(attribute_values.attribute)
            source_values.append(attribute_values.values)
            source_units.append(attribute_values.unit)
            
                
        other.set_values_of_particles_in_units(self.particle_keys, source_attributes, source_values, source_units)
        
    def remove_particles(self, particles):
        indices = self.get_indices_of(particles)
        
        mapping_from_attribute_to_values_and_unit = self.mapping_from_attribute_to_values_and_unit.copy()
        for attribute, attribute_values in mapping_from_attribute_to_values_and_unit.iteritems():
            attribute_values.values = numpy.delete(attribute_values.values,indices)
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
        
    def _remove_particles(self, keys):
        indices = self.get_indices_of(keys)
        
        for attribute, attribute_values in self.mapping_from_attribute_to_values_and_unit.iteritems():
            attribute_values.values = numpy.delete(attribute_values.values,indices)
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
        
    def reindex(self):
        d = {}
        index = 0
        for particle_key in self.particle_keys:
            d[particle_key] = index
            index += 1
          
        self.mapping_from_particle_to_index = d

    def attributes(self):
        return set(self.mapping_from_attribute_to_values_and_unit.keys())
    
    def _state_attributes(self):
        return self.attributes()
        
    def set_model_time(self, value): 
        model_times = self._convert_model_times(value, len(self.particle_keys))
        for attribute_values in self.mapping_from_attribute_to_values_and_unit.values():
            attribute_values.model_times = model_times
