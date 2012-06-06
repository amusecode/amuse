import numpy
from numpy import ma

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
            storage = InMemoryAttribute.new_attribute(attribute, len(keys), quantity)
            self.mapping_from_attribute_to_quantities[attribute] = storage
            storage.set_values(None, quantity)
         
        self.particle_keys = numpy.array(keys, dtype='uint64')

        self.reindex()
    
    def append_to_storage(self, keys, attributes, values):
        for attribute, values_to_set in zip(attributes, values):
            if attribute in self.mapping_from_attribute_to_quantities:
                storage = self.mapping_from_attribute_to_quantities[attribute]
            else:
                storage = InMemoryAttribute.new_attribute(attribute, len(self.particle_keys), values_to_set)
                self.mapping_from_attribute_to_quantities[attribute] = storage
        
            storage.increase_to_length(len(self.particle_keys) + len(keys))
            storage.set_values(slice(len(self.particle_keys), len(self.particle_keys) + len(keys)), values_to_set)
            
        
        old_length = len(self.particle_keys)
        for attribute, attribute_values in list(self.mapping_from_attribute_to_quantities.iteritems()):
            attribute_values.increase_to_length(len(self.particle_keys) + len(keys))
                
        self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys), dtype='uint64')))
        self.reindex()
            
    def get_values_in_store(self, particles, attributes, by_key = True):
        if by_key:
            indices = self.get_indices_of(particles)
        else:
            indices = particles
            
        results = []
        for attribute in attributes:
            storage = self.mapping_from_attribute_to_quantities[attribute]
        
            selected_values = storage.get_values(indices)
        
            results.append(selected_values)
        
        return results
        
    def set_values_in_store(self, particles, attributes, list_of_values_to_set, by_key = True):
        if by_key:
            indices = self.get_indices_of(particles)
        else:
            indices = particles
        
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
    
            if attribute in self.mapping_from_attribute_to_quantities:
                storage = self.mapping_from_attribute_to_quantities[attribute]
            else:
                storage = InMemoryAttribute.new_attribute(attribute, len(self.particle_keys), values_to_set)
                self.mapping_from_attribute_to_quantities[attribute] = storage
            
            try:
                storage.set_values(indices, values_to_set)
            except ValueError as ex:
                # hack to set values between 
                # with quanities with units.none
                # and when values are stored without units
                # need to be removed when units.none is completely gone
                if is_quantity(values_to_set) and not storage.has_units():
                    if not values_to_set.unit.base:
                        storage.set_values(indices, values_to_set.value_in(units.none))
                    else:
                        raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex)) 
                elif not is_quantity(values_to_set) and storage.has_units():
                    if not storage.quantity.unit.base:
                        storage.set_values(indices, units.none.new_quantity(values_to_set))
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
        
        return attribute_values.get_value(index)
   
    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
        
        indices = numpy.searchsorted(self.sorted_keys, particles)
        return self.sorted_indices[indices]
        
    def remove_particles_from_store(self, keys):
        indices = self.get_indices_of(keys)
        
        for attribute, attribute_values in list(self.mapping_from_attribute_to_quantities.iteritems()):
            attribute_values.remove_indices(indices)
            
        
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
        
        return attribute_values.get_value(index)
    
    
class InMemoryGridAttributeStorage(object):
    
    def __init__(self, *number_of_points_in_each_direction):
        self.mapping_from_attribute_to_quantities = {}
        self.number_of_points_in_each_direction = number_of_points_in_each_direction
        
    def can_extend_attributes(self):
        return True
        
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
        copy = InMemoryGridAttributeStorage(*self.number_of_points_in_each_direction)
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
        
        return attribute_values.get_value(index)
   
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
            index = particle_key
        
        return attribute_values.get_value(index)


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
    def copy(self):
        pass

    
    @classmethod
    def _determine_shape(cls, length, input):
        vector_shape = input.shape
        if len(vector_shape) > 1 and len(input) == length:
            return vector_shape
        else:
            return length
            
    @classmethod
    def new_attribute(cls, name, shape, input):
        if is_quantity(input):
            if input.is_vector() :
                shape = cls._determine_shape(shape, input)
            return InMemoryVectorQuantityAttribute(name, shape, input.unit)
        elif hasattr(input, 'as_set'):
            return InMemoryLinkedAttribute(name, shape, input.as_set()._original_set())
        else:
            array = numpy.asanyarray(input)
            dtype = array.dtype
            shape = cls._determine_shape(shape, array)
            return InMemoryUnitlessAttribute(name, shape, dtype)

    def get_value(self, index):
        pass

    def remove_indices(self, indices):
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
        try:
            self.quantity[indices] = values
        except AttributeError:
            if not is_quantity(values):
                raise ValueError("Tried to set a non quantity value for an attribute ({0}) with a unit".format(self.name))
            else:
                raise
    
    def get_shape(self):
        return self.quantity.shape
    

    def increase_to_length(self, newlength):
        delta = newlength - len(self.quantity)
        if delta == 0: 
           return
        deltashape = list(self.quantity.shape)
        deltashape[0] = delta
    
        zeros_for_concatenation = VectorQuantity.zeros(deltashape, self.quantity.unit)
        self.quantity.extend(zeros_for_concatenation)

    def get_length(self):
        return len(self.quantity)

    def copy(self):
        result = type(self)(self.name, self.get_shape(), self.quantity.unit)
        result.set_values(None, self.get_values(None))
        return result

    def get_value(self, index):
        return self.quantity[index]

    def remove_indices(self, indices):
        self.quantity._number = numpy.delete(self.quantity.number, indices)

    def has_units(self):
        return True

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
        if delta == 0: 
           return
        deltashape = list(self.values.shape)
        deltashape[0] = delta
        zeros_for_concatenation =  numpy.zeros(deltashape, dtype = self.values.dtype)
        self.values = numpy.concatenate([self.values, zeros_for_concatenation])

    def get_shape(self):
        return self.values.shape

    def copy(self):
        result = type(self)(self.name, self.get_shape(), self.values.dtype)
        result.set_values(None, self.get_values(None))
        return result

    def get_value(self, index):
        return self.values[index]

    def remove_indices(self, indices):
        self.values = numpy.delete(self.values, indices)

    def has_units(self):
        return False


class InMemoryLinkedAttribute(InMemoryAttribute):
    
    def __init__(self, name, shape, linked_set):
        InMemoryAttribute.__init__(self, name)
        
        self.linked_set = linked_set
        
        self.values = numpy.ma.masked_all(
            shape,
            dtype = linked_set.get_all_keys_in_store().dtype
        )
        
        
    def get_values(self, indices):
        return self.linked_set._masked_subset(self.values[indices])
    
    def set_values(self, indices, values):
        if hasattr(values, 'get_all_keys_in_store'):
            keys = values.get_all_keys_in_store()
            mask = ~values.get_valid_particles_mask()
            keys = numpy.ma.array(keys, dtype=self.values.dtype)
            keys.mask = mask
            self.values[indices] = keys
        else:
            if values is None:
                self.values[indices] = ma.masked
            else:
                for index, key in zip(indices, values):
                    if key is None:
                        self.values[index] = ma.masked
                    else:
                        self.values[index] = key
    
    def get_length(self):
        return self.values.shape
        
    def increase_to_length(self, newlength):
        delta = newlength - len(self.values)
        if delta == 0: 
           return
        deltashape = list(self.values.shape)
        deltashape[0] = delta
        zeros_for_concatenation =  numpy.ma.masked_all(deltashape, dtype = self.values.dtype)
        self.values = numpy.ma.concatenate([self.values, zeros_for_concatenation])

    def get_shape(self):
        return self.values.shape

    def copy(self):
        result = type(self)(self.name, self.get_shape(), self.linked_set)
        result.set_values(None, self.get_values(None))
        return result

    def get_value(self, index):
        key = self.values[index]
        if key is ma.masked:
            return None
        else:
            return self.linked_set._get_particle_unsave(key)

    def remove_indices(self, indices):
        newmask = numpy.delete(self.values.mask, indices)
        newvalues = numpy.delete(self.values.data, indices)
        
        self.values = numpy.ma.array(newvalues, mask = newmask) 

    def has_units(self):
        return False


