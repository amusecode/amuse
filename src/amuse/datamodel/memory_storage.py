import numpy

from amuse.units import quantities
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import is_quantity
from amuse.datamodel.base import *
from amuse.support import exceptions

class InMemoryAttributeStorage(AttributeStorage):
    
    def __init__(self):
        self.mapping_from_attribute_to_quantities = {}
        self.particle_keys = numpy.zeros(0)
        self.__version__ = 0
        
        self.sorted_keys = numpy.zeros(0, dtype=numpy.int32)
        self.sorted_indices = numpy.zeros(0)
        self.index_array = numpy.zeros(0, dtype=numpy.int32)
        self.keys_set = set([])
    
    def can_extend_attributes(self):
        return True
        
    def add_particles_to_store(self, keys, attributes = [], quantities = []):
        if len(quantities) != len(attributes):
            raise exceptions.AmuseException(
                "you need to provide the same number of quantities as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(quantities)
                )
            )
        if len(quantities) > 0 and len(keys) != len(quantities[0]):
            raise exceptions.AmuseException(
                "you need to provide the same number of values as particles, found {0} values and {1} particles".format(
                    len(quantities[0]), len(keys)
                )
            )
        
        if len(self.particle_keys) > 0:
            self.append_to_storage(keys, attributes, quantities)
        else:
            self.setup_storage(keys, attributes, quantities)
    
        self.__version__ = self.__version__ + 1
        
        self.index_array = numpy.arange(len(self))
            
    def setup_storage(self, keys, attributes, quantities):
        self.mapping_from_attribute_to_quantities = {}
        for attribute, quantity in zip(attributes, quantities):
            if not is_quantity(quantity):
                self.mapping_from_attribute_to_quantities[attribute] = quantity
            else:
                self.mapping_from_attribute_to_quantities[attribute] = quantity.as_vector_quantity()
         
        self.particle_keys = numpy.array(keys, dtype='uint64')

        self.reindex()
    
    def append_to_storage(self, keys, attributes, values):
        for attribute, values_to_set in zip(attributes, values):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                if not is_quantity(values_to_set):
                    attribute_values = numpy.zeros(len(self.particle_keys))
                else:
                    attribute_values = VectorQuantity.zeros(
                        len(self.particle_keys),
                        values_to_set.unit,
                    )
            
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            
            
            if not is_quantity(attribute_values):
                attribute_values = numpy.concat(attribute_values,values_to_set)
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            else:
                attribute_values.extend(values_to_set)
        
        old_length = len(self.particle_keys)
        for attribute_values in self.mapping_from_attribute_to_quantities.values():
            if len(attribute_values) == old_length:
                zeros_for_concatenation = VectorQuantity.zeros(len(keys), attribute_values.unit)
                attribute_values.extend(zeros_for_concatenation)
                
        self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys), dtype='uint64')))
        self.reindex()
            
    def get_values_in_store(self, particles, attributes, by_key = True):
        if by_key:
            indices = self.get_indices_of(particles)
        else:
            indices = particles
            
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            if indices is None:
                selected_values = attribute_values
            else:
                selected_values = attribute_values.take(indices)
            
            results.append(selected_values)
        
        return results
        
    def set_values_in_store(self, particles, attributes, list_of_values_to_set, by_key = True):
        if by_key:
            indices = self.get_indices_of(particles)
        else:
            indices = particles
        
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                if is_quantity(values_to_set):
                    attribute_values = VectorQuantity.zeros(
                        len(self.particle_keys),
                        values_to_set.unit,
                    )
                else:
                    dtype = numpy.asanyarray(values_to_set).dtype
                    attribute_values = numpy.zeros(len(self.particle_keys), dtype = dtype)
                    
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            
            try:
                attribute_values.put(indices, values_to_set)
            except ValueError as ex:
                # hack to set values between 
                # with quanities with units.none
                # and when values are stored without units
                # need to be removed when units.none is completely gone
                if is_quantity(values_to_set) and not is_quantity(attribute_values):
                    if not values_to_set.unit.base:
                        attribute_values.put(indices, values_to_set.value_in(units.none))
                    else:
                        raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex)) 
                elif not is_quantity(values_to_set) and is_quantity(attribute_values):
                    if not attribute_values.unit.base:
                        attribute_values.put(indices, units.none.new_quantity(values_to_set))
                    else:
                        raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex)) 
                else:
                    raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex))
                
    def has_key_in_store(self, key):
        return key in self.keys_set
        
    def get_all_keys_in_store(self):
        return self.particle_keys
        
    def get_all_indices_in_store(self):
        return self.index_array
        
    def __len__(self):
        return len(self.particle_keys)
        
    def copy(self):
        copy = get_in_memory_attribute_storage_factory()()
        copy.sorted_keys = self.sorted_keys.copy()
        copy.sorted_indices = self.sorted_indices.copy()
        copy.keys_set = self.keys_set.copy()
        copy.particle_keys = self.particle_keys.copy()
        copy.index_array = self.index_array.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle_key, attribute):
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
        
        index = self.get_indices_of(particle_key)
        
        return attribute_values[index]
   
    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
        
        
        indices = numpy.searchsorted(self.sorted_keys, particles)
        return self.sorted_indices[indices]
        
    def remove_particles_from_store(self, keys):
        indices = self.get_indices_of(keys)
        
        for attribute, attribute_values in list(self.mapping_from_attribute_to_quantities.iteritems()):
            if is_quantity(attribute_values):
                attribute_values._number = numpy.delete(attribute_values.number, indices)
            else:
                attribute_values = numpy.delete(attribute_values, indices)
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
    
        self.__version__ = self.__version__ + 1
        self.index_array = numpy.arange(len(self))
        
    def reindex(self):
        self.sorted_indices = numpy.argsort(self.particle_keys, kind='mergesort')
        self.sorted_keys = self.particle_keys[self.sorted_indices]       
        self.keys_set = set(self.particle_keys)
        self.index_array = numpy.arange(len(self))
        

    def get_defined_attribute_names(self):
        return sorted(self.mapping_from_attribute_to_quantities.keys())
        
        

    def _get_values_for_indices(self, indices, attributes):
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            if indices is None:
                selected_values = attribute_values
            else:
                selected_values = attribute_values.take(indices)
            
            results.append(selected_values)
        
        return results
    
    

    def get_value_in_store(self, particle_key, attribute, by_key = True):
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
        if by_key:
            indices = self.get_indices_of(particle_key)
        else:
            indices = particle_key
        return attribute_values[index]
    
    
class InMemoryGridAttributeStorage(object):
    
    def __init__(self, *number_of_points_in_each_direction):
        self.mapping_from_attribute_to_quantities = {}
        self.number_of_points_in_each_direction = number_of_points_in_each_direction
        
    def storage_shape(self):
        return self.number_of_points_in_each_direction
        
    def add_particles_to_store(self, keys, attributes = [], quantities = []):
        raise exceptions.AmuseException("adding points to the grid is not implemented")
            
    def remove_particles_from_store(self, keys):
        raise exceptions.AmuseException("removing points from the grid is not implemented")
        
    def get_values_in_store(self, indices, attributes):
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            if indices is None:
                selected_values = attribute_values
            else:
                selected_values = attribute_values[indices]
            
            results.append(selected_values)
        
        return results
        
    def set_values_in_store(self,  indices, attributes, list_of_values_to_set):
        
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                if is_quantity(values_to_set):
                    attribute_values = VectorQuantity.zeros(
                       (self.storage_shape()),
                       values_to_set.unit,
                    )
                else: 
                    dtype = numpy.asanyarray(values_to_set).dtype
                    attribute_values = numpy.zeros((self.storage_shape()), dtype = dtype)
                    
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
                 
            attribute_values[indices] = values_to_set
     
    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_to_index
        
    def get_all_keys_in_store(self):
        return None #numpy.s_[0:self.number_of_i], numpy.s_[0:self.number_of_j], numpy.s_[0:self.number_of_k]
        
    def __len__(self):
        return self.number_of_i * self.number_of_j * self.number_of_k
        
    def copy(self):
        copy = InMemoryGridAttributeStorage()
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
    
    def get_defined_attribute_names(self):
        return sorted(self.mapping_from_attribute_to_quantities.keys())
        
    def _get_writeable_attribute_names(self):
        return self.get_defined_attribute_names()




class InMemoryAttributeStorageUseDictionaryForKeySet(InMemoryAttributeStorage):

    

    def __init__(self):
        InMemoryAttributeStorage.__init__(self)
        self.mapping_from_particle_to_index = {}

    

    def append_to_storage(self, keys, attributes, values):
        for attribute, values_to_set in zip(attributes, values):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                if not is_quantity(values_to_set):
                    attribute_values = numpy.zeros(len(self.particle_keys))
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

            

    

    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_to_index

    def copy(self):
        copy = type(self)()
        copy.mapping_from_particle_to_index = self.mapping_from_particle_to_index.copy()
        copy.particle_keys = self.particle_keys.copy()
        copy.index_array = self.index_array.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy

    def get_value_of(self, particle_key, attribute):
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
        index = self.mapping_from_particle_to_index[particle_key]
        return attribute_values[index]

    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
    
        mapping_from_particle_to_index = self.mapping_from_particle_to_index
        result = numpy.zeros(len(particles),dtype='int32')
        index = 0
        for index, particle_key in enumerate(particles):
            result[index] = mapping_from_particle_to_index[particle_key]
            index += 1
        
        return result

    def reindex(self):
        new_index = {}
        index = 0
        for particle_key in self.particle_keys:
            new_index[particle_key] = index
            index += 1
        
        self.mapping_from_particle_to_index = new_index


class InMemoryAttributeStorageUseSortedKeys(InMemoryAttributeStorage):
    
    def __init__(self):
        InMemoryAttributeStorage.__init__(self)
        
        self.sorted_keys = []
        self.sorted_indices = []
        self.keys_set = set([])
    
    def append_to_storage(self, keys, attributes, values):
        for attribute, values_to_set in zip(attributes, values):
            if attribute in self.mapping_from_attribute_to_quantities:
                attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            else:
                if not is_quantity(values_to_set):
                    attribute_values = numpy.zeros(len(self.particle_keys))
                else:
                    attribute_values = VectorQuantity.zeros(
                        len(self.particle_keys),
                        values_to_set.unit,
                    )
                attribute_values = VectorQuantity.zeros(
                    len(self.particle_keys),
                    values_to_set.unit,
                )
            
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            
            if not is_quantity(attribute_values):
                attribute_values = numpy.concatenate((attribute_values,values_to_set))
                self.mapping_from_attribute_to_quantities[attribute] = attribute_values
            else:
                attribute_values.extend(values_to_set)
            
        
        old_length = len(self.particle_keys)
        for attribute, attribute_values in list(self.mapping_from_attribute_to_quantities.iteritems()):
            if len(attribute_values) == old_length:
                
                if not is_quantity(attribute_values):
                    zeros_for_concatenation = numpy.zeros(len(keys), dtype=attribute_values.dtype)
                    attribute_values = numpy.concatenate((attribute_values,zeros_for_concatenation))
                    self.mapping_from_attribute_to_quantities[attribute] = attribute_values
                else:
                    zeros_for_concatenation = VectorQuantity.zeros(len(keys), attribute_values.unit)
                    attribute_values.extend(zeros_for_concatenation)
                
        self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys), dtype='uint64')))
        self.reindex()
            
    
    def has_key_in_store(self, key):
        return key in self.keys_set
        
    def copy(self):
        copy = type(self)()
        copy.sorted_keys = self.sorted_keys.copy()
        copy.sorted_indices = self.sorted_indices.copy()
        copy.keys_set = self.keys_set.copy()
        copy.particle_keys = self.particle_keys.copy()
        copy.index_array = self.index_array.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle_key, attribute):
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
        
        index = self.get_indices_of(particle_key)
        
        return attribute_values[index]
   
    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
        
        if len(self.particle_keys) == 0:
            return ()
            
        indices = numpy.searchsorted(self.sorted_keys, particles)
        return self.sorted_indices[indices]
        
        
    def reindex(self):
        self.sorted_indices = numpy.argsort(self.particle_keys, kind='mergesort')
        self.sorted_keys = self.particle_keys[self.sorted_indices]       
        self.keys_set = set(self.particle_keys)
        

        
    def get_value_in_store(self, particle_key, attribute, by_key = True):
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
        if by_key:
            index = self.get_indices_of(particle_key)
        else:
            index = particles
        return attribute_values[index]


def get_in_memory_attribute_storage_factory():
    if False:
        return InMemoryAttributeStorageUseDictionaryForKeySet
    else:
        return InMemoryAttributeStorageUseSortedKeys


class InMemoryAttribute(object):
    
    def __init__(self, name):
        self.name = name
        
    def get_values(self, indices):
        pass
    
    def set_values(self, indices, valus):
        pass
    
    def get_length(self):
        return 0
    
    def get_shape(self):
        return 0
    
    def increase_to_length(self, newlength):
        pass
class InMemoryVectorQuantityAttribute(InMemoryAttribute):
    
    def __init__(self, name, shape, unit):
        InMemoryAttribute.__init__(self, name)
    
        self.quantity = VectorQuantity.zeros(
            shape,
            unit,
        )
        
    def get_values(self, indices):
        return self.quantity[indices]
    
    def set_values(self, indices, values):
        self.quantity[indices] = values
    
    def get_shape(self):
        return self.quantity.shape
    

    def increase_to_length(self, newlength):
        delta = newlength - len(self.quantity)
        deltashape = list(self.quantity.shape)
        deltashape[0] = delta
    
        zeros_for_concatenation = VectorQuantity.zeros(deltashape, self.quantity.unit)
        self.quantity.extend(zeros_for_concatenation)

    def get_length(self):
        return len(self.quantity)

class InMemoryUnitlessAttribute(InMemoryAttribute):
    
    def __init__(self, name, shape, dtype = 'float64'):
        InMemoryAttribute.__init__(self, name)
        
        self.values = numpy.zeros(
            shape,
            dtype = dtype
        )
        
    def get_values(self, indices):
        return self.values[indices]
    
    def set_values(self, indices, values):
        self.values[indices] = values
    
    def get_length(self):
        return self.values.shape
    

    def increase_to_length(self, newlength):
        delta = newlength - len(self.values)
        deltashape = list(self.values.shape)
        deltashape[0] = delta

        zeros_for_concatenation =  numpy.zeros(deltashape, dtype = self.values.dtype)
        self.values = numpy.concatenate(self.values ,zeros_for_concatenation)

    def get_shape(self):
        return self.quantity.shape

