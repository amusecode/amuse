import numpy
from numpy import ma

from itertools import izip

from amuse.units import quantities
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import is_quantity
from amuse.datamodel.base import *
from amuse.support import exceptions
from amuse.rfi.async_request import FakeASyncRequest

try:
  if numpy.uintp().itemsize<8: raise Exception()
  from amuse.datamodel.simple_hash import SimpleHash
  _SIMPLE_HASH_PRESENT_=True
except BaseException as ex:
  _SIMPLE_HASH_PRESENT_=False

_PREFER_SORTED_KEYS_=True

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
        
    
    def remove_attribute_from_store(self, name):
        del self.mapping_from_attribute_to_quantities[name]
        
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
        
        self.__version__ = self.__version__ + 1
        self.index_array = numpy.arange(len(self.particle_keys) + len(keys))
        
        if len(self.particle_keys) > 0:
            previous_length = len(self.particle_keys)
            self.append_to_storage(keys, attributes, quantities)
            return self.index_array[previous_length:]
        else:
            self.setup_storage(keys, attributes, quantities)
            return self.index_array
    
        
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
            try:
                storage.set_values(slice(len(self.particle_keys), len(self.particle_keys) + len(keys)), values_to_set)
            except Exception as ex:
                raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex))
                
        old_length = len(self.particle_keys)
        for attribute, attribute_values in list(self.mapping_from_attribute_to_quantities.iteritems()):
            attribute_values.increase_to_length(len(self.particle_keys) + len(keys))
                
        self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys), dtype='uint64')))
        self.reindex()
            
    def get_values_in_store(self, indices, attributes):            
        results = []
        for attribute in attributes:
            if not attribute in self.mapping_from_attribute_to_quantities:
                raise AttributeError('"{0}" not defined for grid'.format(attribute))
                
            storage = self.mapping_from_attribute_to_quantities[attribute]
        
            selected_values = storage.get_values(indices)
        
            results.append(selected_values)
        
        return results
    
    def get_values_in_store_async(self, indices, attributes):
        result = self.get_values_in_store(indices, attributes)
        return FakeASyncRequest(result)
        
    def set_values_in_store(self, indices, attributes, list_of_values_to_set):
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
                # note an alternative might be to store always with units.none
                # but have the getitem remove the unit, caveat: maintaining compatibility in fileformat?
                if is_quantity(values_to_set) and not storage.has_units():
                    if not values_to_set.unit.base:
                        storage.set_values(indices, values_to_set.value_in(units.none))
                    else:
                        raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex)) 
                # this is no longer necessary:
                #~ elif not is_quantity(values_to_set) and storage.has_units():
                    #~ if not storage.quantity.unit.base:
                        #~ storage.set_values(indices, units.none.new_quantity(values_to_set))
                    #~ else:
                        #~ raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex)) 
                else:
                    raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex))

    def set_values_in_store_async(self, indices, attributes, list_of_values_to_set):
        result = self.set_values_in_store(indices, attributes, list_of_values_to_set)
        return FakeASyncRequest(result)
                
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
        
    def get_value_of(self, index, attribute):
        return self.get_value_in_store(index, attribute)
   
    def get_indices_of(self, keys):
        raise NotImplementedError()
        
    def remove_particles_from_store(self, indices):
            
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
        
    def get_defined_settable_attribute_names(self):
        return self.get_defined_attribute_names()
        

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


    def get_value_in_store(self, index, attribute):
        attribute_values = self.mapping_from_attribute_to_quantities[attribute]
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
        
    def remove_attribute_from_store(self, name):
        del self.mapping_from_attribute_to_quantities[name]
        
        
    def get_values_in_store(self, indices, attributes):
        results = []
        for attribute in attributes:
            if not attribute in self.mapping_from_attribute_to_quantities:
                raise AttributeError('"{0}" not defined for grid'.format(attribute))
            attribute_values = self.mapping_from_attribute_to_quantities[attribute]
            #if indices is None:
            #    selected_values = attribute_values
            #else:
            #    selected_values = attribute_values[indices]
            
            results.append(attribute_values.get_values(indices))
        
        return results
        
    def set_values_in_store(self,  indices, attributes, list_of_values_to_set):
        
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):

            if attribute in self.mapping_from_attribute_to_quantities:
                storage = self.mapping_from_attribute_to_quantities[attribute]
            else:
                storage = InMemoryAttribute.new_attribute(attribute, self.storage_shape(), values_to_set)
                self.mapping_from_attribute_to_quantities[attribute] = storage
                #storage.set_values(None, quantity)
                #if is_quantity(values_to_set):
                #    attribute_values = VectorQuantity.zeros(
                #       (self.storage_shape()),
                #       values_to_set.unit,
                #    )
                #else: 
                #    dtype = numpy.asanyarray(values_to_set).dtype
                #    attribute_values = numpy.zeros((self.storage_shape()), dtype = dtype)
                    
            try:
                storage.set_values(indices, values_to_set)
            except ValueError as ex:
                # hack to set values between 
                # with quanities with units.none
                # and when values are stored without units
                # note an alternative might be to store always with units.none
                # but have the getitem remove the unit, caveat: maintaining compatibility in fileformat?
                if is_quantity(values_to_set) and not storage.has_units():
                    if not values_to_set.unit.base:
                        storage.set_values(indices, values_to_set.value_in(units.none))
                    else:
                        raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex)) 
                else:
                    raise AttributeError("exception in setting attribute '{0}', error was '{1}'".format(attribute, ex))

     
    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_to_index
        
    def get_all_keys_in_store(self):
        return Ellipsis #numpy.s_[0:self.number_of_i], numpy.s_[0:self.number_of_j], numpy.s_[0:self.number_of_k]
        
    def __len__(self):
        return self.storage_shape()[0]
        
    def copy(self):
        copy = InMemoryGridAttributeStorage(*self.number_of_points_in_each_direction)
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
    
    def get_defined_attribute_names(self):
        return sorted(self.mapping_from_attribute_to_quantities.keys())
        
    def _get_writeable_attribute_names(self):
        return self.get_defined_attribute_names()
    
    def get_defined_settable_attribute_names(self):
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


    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
    
        mapping_from_particle_to_index = self.mapping_from_particle_to_index
        result = []
        notfoundkeys = []
        foundkeys = []
        for index, particle_key in enumerate(particles):
            try:
                result.append(mapping_from_particle_to_index[particle_key])
                foundkeys.append(particle_key)
            except KeyError:
                notfoundkeys.append(particle_key)
        
        if not len(notfoundkeys) == 0:
            raise exceptions.KeysNotInStorageException(
                numpy.asarray(foundkeys), 
                numpy.asarray(result), 
                numpy.asarray(notfoundkeys)
            )
        return numpy.asarray(result)

    def reindex(self):
        new_index=dict(izip(self.particle_keys,xrange(len(self.particle_keys))))
        self.mapping_from_particle_to_index = new_index


class InMemoryAttributeStorageUseSortedKeys(InMemoryAttributeStorage):
    
    def __init__(self):
        InMemoryAttributeStorage.__init__(self)
        
        self.sorted_keys = []
        self.sorted_indices = []
    
    def has_key_in_store(self, key):
        i=numpy.searchsorted(self.sorted_keys, key)
        return i<len(self.sorted_keys) and self.sorted_keys[i]==key
        
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
   
    def get_indices_of(self, keys):
        if keys is None:
            return self.index_array

        if len(self.particle_keys) == 0:
            return ()
        
        indices = numpy.searchsorted(self.sorted_keys, keys)  
        indices = numpy.where(indices >= len(self.sorted_keys), 0, indices)  
        foundkeys = self.sorted_keys[indices]
        are_found = foundkeys == keys
        are_all_found = numpy.all(are_found)
        if not are_all_found:
            arrayedkeys = numpy.asarray(keys)
            notfoundkeys = arrayedkeys[numpy.logical_not(are_found)]
            raise exceptions.KeysNotInStorageException(
                arrayedkeys[are_found], 
                self.sorted_indices[indices[are_found]], 
                notfoundkeys
            )
        return self.sorted_indices[indices]
        
        
    def reindex(self):
        self.sorted_indices = numpy.argsort(self.particle_keys, kind='mergesort')
        self.sorted_keys = self.particle_keys[self.sorted_indices]
        
class InMemoryAttributeStorageUseSimpleHash(InMemoryAttributeStorage):
    
    def __init__(self):
        InMemoryAttributeStorage.__init__(self)
        self._hash=SimpleHash()
    
    def has_key_in_store(self, key):
        return self._hash.key_present(key)
        
    def copy(self):
        copy = type(self)()
        copy.particle_keys = self.particle_keys.copy()
        copy.index_array = self.index_array.copy()
        copy._hash.reindex(copy.particle_keys)
        for attribute, attribute_values in self.mapping_from_attribute_to_quantities.iteritems():
            copy.mapping_from_attribute_to_quantities[attribute] = attribute_values.copy()
        return copy
   
    def get_indices_of(self, keys):
        if keys is None:
            return self.index_array

        if len(self.particle_keys) == 0:
            return ()
        
        return self._hash.lookup(keys)
        
    def reindex(self):
        self._hash.reindex(self.particle_keys)

    def __getstate__(self):
        state=self.__dict__.copy()
        state.pop("_hash")
        return state
      

    def __setstate__(self,state):
        self.__dict__ = state
        self._hash=SimpleHash()
        if len(self.particle_keys):
          self._hash.reindex(self.particle_keys)



def get_in_memory_attribute_storage_factory():
    if _SIMPLE_HASH_PRESENT_:
       return InMemoryAttributeStorageUseSimpleHash
    elif _PREFER_SORTED_KEYS_:
        return InMemoryAttributeStorageUseSortedKeys
    else:   
        return InMemoryAttributeStorageUseDictionaryForKeySet



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
    def _determine_shape(cls, length, values_to_set):
        if isinstance(length, tuple):
            vector_shape = values_to_set.shape
            if len(vector_shape) > len(length):
                return vector_shape
            else:
                return length
        vector_shape = values_to_set.shape
        if len(vector_shape) > 1 and len(values_to_set) == length:
            return vector_shape
        else:
            return length

            
    @classmethod
    def new_attribute(cls, name, shape, values_to_set):
        if is_quantity(values_to_set):
            if values_to_set.is_vector() :
                shape = cls._determine_shape(shape, values_to_set)
            return InMemoryVectorQuantityAttribute(name, shape, values_to_set.unit)
        elif values_to_set is None:
            return InMemoryLinkedAttribute(name, shape)
        else:
            array = numpy.asanyarray(values_to_set)
            dtype = array.dtype
            shape = cls._determine_shape(shape, array)
            if dtype.kind == 'S' or dtype.kind == 'U':
                return InMemoryStringAttribute(name, shape, dtype)
            elif dtype == numpy.object:
                return InMemoryLinkedAttribute(name, shape)
            else:
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
        if indices is None:
            return self.quantity
        else:
            return self.quantity[indices]
    

    def set_values(self, indices, values):
        try:
            if indices is None:
                indices = slice(None)
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
        if indices is None:
            return self.values
        else:
            return self.values[indices]
    
    def set_values(self, indices, values):
        if indices is None:
            indices = slice(None)
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


class InMemoryStringAttribute(InMemoryUnitlessAttribute):
    
    def set_values(self, indices, values):
        if indices is None:
            indices = slice(None)
        if isinstance(values, basestring):
            dtype=numpy.dtype(self.values.dtype.kind+str(max(1,len(values))))
        else:
            values_as_array = numpy.asarray(values, dtype=numpy.dtype(self.values.dtype.kind))
            dtype = values_as_array.dtype
       
        if dtype.itemsize > self.values.dtype.itemsize:
            self.values = numpy.asarray(self.values, dtype=dtype)
            self.dtype = dtype
                        
        self.values[indices] = values
    


        

class InMemoryLinkedAttribute(InMemoryAttribute):
    
    def __init__(self, name, shape):
        InMemoryAttribute.__init__(self, name)
        self.values = LinkedArray(numpy.empty(
            shape,
            dtype = numpy.object
        ))
        
        
    def get_values(self, indices):
        if indices is None:
            objects = self.values
        else:
            objects = self.values[indices]
            
        return objects
                        
    def set_values(self, indices, values):
        if indices is None:
            indices = slice(None)
        self.values[indices] = values
            

    def get_length(self):
        return self.values.shape
        
    def increase_to_length(self, newlength):
        delta = newlength - len(self.values)
        if delta == 0: 
            return
        deltashape = list(self.values.shape)
        deltashape[0] = delta
        zeros_for_concatenation =  numpy.empty(deltashape, dtype = self.values.dtype)
        self.values = LinkedArray(numpy.concatenate([self.values, zeros_for_concatenation]))



    def get_shape(self):
        return self.values.shape

    def copy(self):
        result = type(self)(self.name, self.get_shape())
        result.set_values(None, self.get_values(None))
        return result

    def get_value(self, index):
        value = self.values[index]
        return value

    def remove_indices(self, indices):
        self.values = numpy.delete(self.values, indices)

    def has_units(self):
        return False

