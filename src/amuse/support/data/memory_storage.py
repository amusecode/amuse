import numpy

from amuse.support.data import values
from amuse.support.data.values import VectorQuantity
from amuse.support.data.base import *
    

class InMemoryAttributeStorage(AttributeStorage):
    
    def __init__(self):
        self.mapping_from_attribute_to_quantities = {}
        self.mapping_from_particle_to_index = {}
        self.particle_keys = []

    def _add_particles(self, keys, attributes = [], quantities = []):
        if len(quantities) != len(attributes):
            raise Exception(
                "you need to provide the same number of quantities as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(quantities)
                )
            )
        if len(quantities) > 0 and len(keys) != len(quantities[0]):
            raise Exception(
                "you need to provide the same number of values as particles, found {0} values and {1} particles".format(
                    len(quantities[0]), len(keys)
                )
            )
        
        if len(self.particle_keys) > 0:
            self.append_to_storage(keys, attributes, quantities)
        else:
            self.setup_storage(keys, attributes, quantities)
            
    def setup_storage(self, keys, attributes, quantities):
        self.mapping_from_attribute_to_quantities = {}
        for attribute, quantity in zip(attributes, quantities):
            self.mapping_from_attribute_to_quantities[attribute] = quantity.as_vector_quantity()
         
        self.particle_keys = numpy.array(keys, dtype='uint64')

        self.reindex()
    
    def append_to_storage(self, keys, attributes, values):
        for attribute, values_to_set in zip(attributes, values):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                attribute_values = VectorQuantity.zeros(
                    len(self.particle_keys),
                    values_to_set.unit,
                )
            
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            attribute_values.extend(values_to_set)
        
        old_length = len(self.particle_keys)
        for attribute_values in self.mapping_from_attribute_to_quantities.values():
            if len(attribute_values) == old_length:
                zeros_for_concatenation = VectorQuantity.zeros(len(keys), attribute_values.unit)
                attribute_values.extend(zeros_for_concatenation)
        
                
        index = len(self.particle_keys)

        self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys), dtype='uint64')))

        for particle_key in keys:
            self.mapping_from_particle_to_index[particle_key] = index
            index += 1
            
    def _get_values(self, particles, attributes):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            if indices is None:
                selected_values = attribute_values
            else:
                selected_values = attribute_values.take(indices)
            
            results.append(selected_values)
        
        return results
        
    def _set_values(self, particles, attributes, list_of_values_to_set):
        indices = self.get_indices_of(particles)
        
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                attribute_values = VectorQuantity.zeros(
                   len(self.particle_keys),
                   values_to_set.unit,
                )
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
                 
            attribute_values.put(indices, values_to_set)
            
            
    
    
    def _get_attribute_names(self):
        return sorted(self.mapping_from_attribute_to_quantities.keys())
    
    
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
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle_key, attribute):
        if not attribute in self.mapping_from_attribute_to_quantities:
            raise AttributeError("particle does not have a "+attribute)
        
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
        
        index = self.mapping_from_particle_to_index[particle_key]
        
        return attribute_values[index]
        
            
            
   
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
        
    def _remove_particles(self, keys):
        indices = self.get_indices_of(keys)
        
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            attribute_values._number = numpy.delete(attribute_values.number, indices)
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
        
    def reindex(self):
        new_index = {}
        index = 0
        for particle_key in self.particle_keys:
            new_index[particle_key] = index
            index += 1
          
        self.mapping_from_particle_to_index = new_index

    def attributes(self):
        return set(self.mapping_from_attribute_to_quantities.keys())
    
    def _state_attributes(self):
        return self.attributes()
        
        


class InMemoryGridAttributeStorage(object):
    
    def __init__(self, number_of_i, number_of_j, number_of_k):
        self.mapping_from_attribute_to_quantities = {}
        self.number_of_i = number_of_i
        self.number_of_j = number_of_j
        self.number_of_k = number_of_k
        
    def storage_shape(self):
        return (self.number_of_i, self.number_of_j, self.number_of_k)
        
    def _add_particles(self, keys, attributes = [], quantities = []):
        raise Exception("adding points to the grid is not implemented")
            
    def _remove_particles(self, keys):
        raise Exception("removing points from the grid is not implemented")
        
    def _get_values(self, indices, attributes):
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            if indices is None:
                selected_values = attribute_values
            else:
                selected_values = attribute_values[indices]
            
            results.append(selected_values)
        
        return results
        
    def _set_values(self,  indices, attributes, list_of_values_to_set):
        
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                attribute_values = VectorQuantity.zeros(
                   (self.storage_shape()),
                   values_to_set.unit,
                )
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
                 
            attribute_values[indices] = values_to_set
     
    def _get_attribute_names(self):
        return sorted(self.mapping_from_attribute_to_quantities.keys())
        
    def _has_key(self, key):
        return key in self.mapping_from_particle_to_index
        
    def _get_keys(self):
        return None #numpy.s_[0:self.number_of_i], numpy.s_[0:self.number_of_j], numpy.s_[0:self.number_of_k]
        
    def __len__(self):
        return self.number_of_i * self.number_of_j * self.number_of_k
        
    def copy(self):
        copy = InMemoryGridAttributeStorage()
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
        
    def attributes(self):
        return set(self.mapping_from_attribute_to_quantities.keys())
    
    def _state_attributes(self):
        return self.attributes()
        
    def _get_writeable_attribute_names(self):
        return self.attributes()

