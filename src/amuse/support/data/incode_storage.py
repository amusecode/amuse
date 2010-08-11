from amuse.support.data import parameters
from amuse.support.data.core import Particles, ParticleInformationChannel, Particle
from amuse.support.data.core import AttributeStorage
from amuse.support.methods import AbstractCodeMethodWrapper
import numpy

from amuse.support.units import nbody_system, units
from amuse.support.data import values
from amuse.support.data import base
from amuse.support.core import late
from amuse.support import exceptions

import inspect


class ParticleMappingMethod(AbstractCodeMethodWrapper):
    def __init__(self, method, attribute_names = None):
        AbstractCodeMethodWrapper.__init__(self, method)
        
        if attribute_names is None:
            self._attribute_names = []
        else:
            self._attribute_names = attribute_names
        
    @late
    def name_of_the_indexing_parameter(self):
        return 'index_of_the_particle'
        

class ParticleGetAttributesMethod(ParticleMappingMethod):
    """
    Instances wrap other methods and provide mappings
    from attribute names to results.
    
    Simple attribute getter methods take an array of indices
    and return a tuple with arrays of result values.
    
    .. code::
        x, y, z = instance.get_xyz(indices)
    
    Instances of this class make it possible to access the 
    return values by their attribute names.
    
    For this it employs two strategies:
    
    1. It uses the provided array of names and
       maps each name to the positional output.
    
    2. If no array of names is provided it asks the wrapped 
       method for all the names of the output parameters 
       (this scheme only works for legacy 
       functions or for wrapped legacy functions)
    
    """
    def __init__(self, method, attribute_names = None):
        ParticleMappingMethod.__init__(self, method, attribute_names)
    
    @late
    def attribute_names(self):
        if self._attribute_names:
            return self._attribute_names
        else:
            result = []
            for x in self.method_output_argument_names:
                if x == self.name_of_the_indexing_parameter:
                    continue
                else:
                    result.append(x)
            return result
    
    def check_arguments(self, storage, attributes_to_return, *indices):
        if len(indices[0]) > 1: 
            if self.method_is_legacy and not self.method.specification.can_handle_array:
                raise Exception(
                    "getter method {0} cannot handle arrays".format(self.method)
                )
            elif self.method_is_code:
                if not self.method.legacy_specification is None:
                    if not self.method.legacy_specification.can_handle_array:
                        raise exceptions.AmuseException(
                            "getter method {0} cannot handle arrays".format(self.method)
                        )
    
    def convert_return_value(self, return_value, storage, attributes_to_return):
        if len(self.attribute_names) == 1:
            return_value = (return_value,)
        
        set_of_attributes_to_return = set(attributes_to_return)
        
        result = {}
        
        if self.index_output_attributes:
            index_output_attributes = self.index_output_attributes
        else:
            index_output_attributes = [False] * len(return_value)
        
        for value, attribute, isindex in zip(return_value, self.attribute_names, index_output_attributes):
            if attribute in set_of_attributes_to_return:
                if isindex:
                    result[attribute] = values.new_quantity(storage._get_keys_for_indices_in_the_code(value), units.object_key)
                else:
                    result[attribute] = value
                    
        return result
    
    def get_attribute_values(self, storage, attributes_to_return, *indices):
        
        self.check_arguments(storage, indices, attributes_to_return)
        
        return_value = self.method(*indices, **storage.extra_keyword_arguments_for_getters_and_setters)
        
        return self.convert_return_value(return_value, storage, attributes_to_return)
    
class ParticleSetAttributesMethod(ParticleMappingMethod):
    """
    Instances wrap other methods and provide mappings
    from attribute names to input parameters.
    
    Simple attribute setter methods take an array of indices
    and one or more arrays of new values.
    
    .. code::
       instance.set_xyz(indices, x, y, z)
    
    Instances of this class make it possible to access the 
    possitional parameters with attribute names.
    
    .. Note::
        the index argument is assumed to always come first!
    
    For this it employs two strategies:
    
    1. It uses the provided array of names and
       maps each name to the positional output.
    
    2. If no array of names is provided it asks the wrapped 
       method for all the names of the input parameters 
       (this scheme works for legacy 
       functions and sometimes for python native functions (if
       they have named arguments))
    
    """
    def __init__(self, method,  attribute_names = None):
        ParticleMappingMethod.__init__(self, method,  attribute_names)
    
    @late
    def attribute_names(self):
        if self._attribute_names:
            return self._attribute_names
        else:
            result = []
            for x in self.method_input_argument_names:
                if x == self.name_of_the_indexing_parameter:
                    continue
                else:
                    result.append(x)
            return result
    
    @late
    def names_to_index(self):
        result = {}
        for index, name in enumerate(self.attribute_names):
            result[name] = index
        return result
        
    def set_attribute_values(self, storage, attributes, values, *indices):
        list_arguments = list(indices)
        list_arguments.extend(self.convert_attributes_and_values_to_list_arguments(attributes, values))
        self.method(*list_arguments, **storage.extra_keyword_arguments_for_getters_and_setters)
    
    def convert_attributes_and_values_to_list_arguments(self, attributes, values):
        list_arguments = [0] * (len(self.attribute_names))
        
        names_to_index = self.names_to_index
        for attribute, quantity in zip(attributes, values):
            if attribute in names_to_index:
                index = names_to_index[attribute]
                list_arguments[index] = quantity
        
        return list_arguments


class NewParticleMethod(ParticleSetAttributesMethod):
    
    def __init__(self,  method, attribute_names = None):
        ParticleSetAttributesMethod.__init__(self, method, attribute_names)

    def add_entities(self, attributes, values):
        list_arguments = self.convert_attributes_and_values_to_list_arguments(attributes, values)
        indices = self.method(*list_arguments)
        return indices
        
class ParticleQueryMethod(object):
    def __init__(self, method, names = (), public_name = None):
        self.method = method
        self.name_of_the_out_parameters = names
        self.public_name = public_name

    def apply(self, particles, *args, **kwargs):
        indices = self.method(*args, **kwargs)
        
        keys = particles._private.attribute_storage._get_keys_for_indices_in_the_code(indices)
        
        return particles._subset(keys)
        

class ParticleSpecificSelectMethod(object):
    def __init__(self, method, names = (), public_name = None):
        self.method = method
        self.name_of_the_out_parameters = names
        self.public_name = public_name

    def apply_on_all(self, particles):
        
        all_indices = particles._private.attribute_storage.mapping_from_index_in_the_code_to_particle_key.keys()
        
        lists_of_indices = self.method(list(all_indices))
        
        lists_of_keys = []
        for indices in lists_of_indices:
                keys = particles._private.attribute_storage._get_keys_for_indices_in_the_code(indices)        
                lists_of_keys.append(keys)
        
        result = []
        for keys in zip(list_of_keys):
            result.append(particles._subset(keys))
            
        return result
    
    def apply_on_one(self, set,  particle):
        
        index = set._private.attribute_storage.get_indices_of(particle.key)
        
        result = self.method(index)
        
        keys = set._private.attribute_storage._get_keys_for_indices_in_the_code(result)  
        
        result = []
        return particles._subset(keys)
        
        
class ParticleMethod(AbstractCodeMethodWrapper):
    
    def __init__(self, method, public_name = None):
        AbstractCodeMethodWrapper.__init__(self, method)
        self.public_name = public_name

    def apply_on_all(self, particles, *list_arguments, **keyword_arguments):
        storage = particles._private.attribute_storage
        all_indices = storage.mapping_from_index_in_the_code_to_particle_key.keys()
        return self.method(list(all_indices), *list_arguments, **keyword_arguments)
    
    def apply_on_one(self, set,  particle, *list_arguments, **keyword_arguments):
        storage = particle.particles_set._private.attribute_storage
        index = storage.get_indices_of([particle.key])
        return self.method(index[0], *list_arguments, **keyword_arguments)
        



class AbstractInCodeAttributeStorage(base.AttributeStorage):
    def __init__(self, 
            code_interface, 
            setters,
            getters,
            extra_keyword_arguments_for_getters_and_setters = {},
):
        
        self.code_interface = code_interface
        
        self.getters = getters
        self.setters = setters
        
        self.attributes = set([])
        for x in self.getters:
            self.attributes |= set(x.attribute_names)
        for x in self.setters:
            self.attributes |= set(x.attribute_names)
            
        self.writable_attributes = set([])
        for x in self.setters:
            self.writable_attributes |= set(x.attribute_names)
            
        
        self.extra_keyword_arguments_for_getters_and_setters = extra_keyword_arguments_for_getters_and_setters
        
    
    def select_getters_for(self, attributes):
        set_of_attributes = set(attributes)
        
        # first check for an exact match
        result = [getter for getter in self.getters if set(getter.attribute_names) == set_of_attributes]
        if result:
            return result
        
        # sort methods on attribute lengths, longest first
        sorted_getters = sorted(self.getters, key=lambda x : len(x.attribute_names), reverse = True)
        
        # next, select the longest fitting method(s), to minize the number of calls
        for access_method in sorted_getters:
            if set_of_attributes >= set(access_method.attribute_names):
                result.append(access_method)
                set_of_attributes -= set(access_method.attribute_names)
        
        # next, select the sortest method(s), to minimize the extra parameters
        if set_of_attributes:
            for access_method in reversed(sorted_getters):
                if set_of_attributes & set(access_method.attribute_names):
                    result.append(access_method)
                    set_of_attributes -= set(access_method.attribute_names)
                    
        if set_of_attributes:
            raise exceptions.AmuseException("Do not have attributes {0}".format(sorted(set_of_attributes)))
        
        return result
    
    def select_setters_for(self, attributes):
        set_of_attributes = set(attributes)
        result = []
        for access_method in self.setters:
            if set_of_attributes >= set(access_method.attribute_names):
                result.append(access_method)
                set_of_attributes -= set(access_method.attribute_names)
                
        if set_of_attributes:
            raise exceptions.AmuseException("Cannot set attributes {0}".format(sorted(set_of_attributes)))
            
        return result
    
    def _get_attribute_names(self):
        return self.attributes
    
    def _state_attributes(self):
        return self._get_attribute_names()
        
        

class InCodeAttributeStorage(AbstractInCodeAttributeStorage):
       
    def __init__(self, 
            code_interface, 
            new_particle_method, 
            delete_particle_method, 
            number_of_particles_method, 
            setters,
            getters,
            name_of_the_index):
        
        
        for x in getters:
            x.name_of_the_indexing_parameter = name_of_the_index
            
        for x in setters:
            x.name_of_the_indexing_parameter = name_of_the_index
        
        AbstractInCodeAttributeStorage.__init__(self, code_interface, setters, getters)

        self.mapping_from_particle_key_to_index_in_the_code = {}
        self.mapping_from_index_in_the_code_to_particle_key = {}
        self.particle_keys = []
        
        self._get_number_of_particles = number_of_particles_method
        self.delete_particle_method = delete_particle_method
        self.new_particle_method = new_particle_method
        
        
                    
    def __len__(self):
        return self._get_number_of_particles()
            
            
    def _add_particles(self, keys, attributes = [], values = []):
        
        indices = self.new_particle_method.add_entities(attributes, values)
        
        
        if len(self.particle_keys) > 0:
            self.particle_keys = numpy.concatenate((self.particle_keys, numpy.array(list(keys))))
        else:
            self.particle_keys = numpy.array(keys)

        index = 0
        for key in keys:
            self.mapping_from_particle_key_to_index_in_the_code[key] = indices[index]
            self.mapping_from_index_in_the_code_to_particle_key[indices[index]] = key
            index = index + 1
            
        
        
    def get_indices_of(self, keys):
        indices_in_the_code = []
        if keys is None:
            keys = self.particle_keys
            
        for particle_key in keys:
            indices_in_the_code.append(self.mapping_from_particle_key_to_index_in_the_code[particle_key])
            
        return indices_in_the_code
        
   
    def get_key_indices_of(self, keys):
        result = []
        if keys is None:
            keys = self.particle_keys
        
        keys_set = set(keys)
        for index in range(len(self.particle_keys)):
           key = self.particle_keys[index]
           if key in keys_set:
               result.append(index)
          
        return result
         
    def _get_values(self, keys, attributes):
        indices_in_the_code = self.get_indices_of(keys)
        
        if len(indices_in_the_code) == 0:
            return [[] for attribute in attributes]
             
        mapping_from_attribute_to_result = {}
        
        for getter in self.select_getters_for(attributes):
            result = getter.get_attribute_values(self, attributes, indices_in_the_code)
            mapping_from_attribute_to_result.update(result)
            
        results = []
        for attribute in attributes:
            results.append(mapping_from_attribute_to_result[attribute])
        return results
        
    def _set_values(self, keys, attributes, values):
        indices_in_the_code = self.get_indices_of(keys)
        
        if len(indices_in_the_code) == 0:
            return
            
        for setter in self.select_setters_for(attributes):
            setter.set_attribute_values(self, attributes, values, indices_in_the_code)
    
    def _remove_particles(self, keys):
        indices_in_the_code = self.get_indices_of(keys)
        
        self.delete_particle_method(indices_in_the_code)
        
        d = self.mapping_from_particle_key_to_index_in_the_code
        for key in keys:
            del d[key]
        
        for i in indices_in_the_code:
            del self.mapping_from_index_in_the_code_to_particle_key[i]
        
         
        indices_to_delete = self.get_key_indices_of(keys)
        self.particle_keys =  numpy.delete(self.particle_keys, indices_to_delete)
            
        
    def _get_keys(self):
        return self.particle_keys

    def _has_key(self, key):
        return key in self.mapping_from_particle_key_to_index_in_the_code
        
    def _get_keys_for_indices_in_the_code(self, indices):
        result = []
        for i in indices:
            result.append(self.mapping_from_index_in_the_code_to_particle_key.get(i, -1))
        return result
    

class InCodeGridAttributeStorage(AbstractInCodeAttributeStorage):
    
    def __init__(self, 
            code_interface, 
            get_range_method,
            setters,
            getters,
            extra_keyword_arguments_for_getters_and_setters = {},
    ):
        AbstractInCodeAttributeStorage.__init__(self, code_interface, setters, getters, extra_keyword_arguments_for_getters_and_setters)
        self.get_range_method = get_range_method
            
    def storage_shape(self):
        imin, imax, jmin, jmax, kmin, kmax = self.get_range_method()
        return (imax - imin + 1, jmax - jmin + 1, kmax - kmin + 1)
        
    def _add_particles(self, keys, attributes = [], quantities = []):
        raise exceptions.AmuseException("adding points to the grid is not implemented")
            
    def _remove_particles(self, keys):
        raise exceptions.AmuseException("removing points from the grid is not implemented")
    
    def _to_arrays_of_indices(self, index):
        imin, imax, jmin, jmax, kmin, kmax = self.get_range_method()
        indices = numpy.mgrid[slice(imin, imax+1),slice(jmin, jmax+1),slice(kmin, kmax+1)]
        indices = [x[index] for x in indices]
        return indices
        
    def _get_values(self, indices, attributes):
        i,j,k = self._to_arrays_of_indices(indices)
        mapping_from_attribute_to_result = {}
        for getter in self.select_getters_for(attributes):
            result = getter.get_attribute_values(self, attributes, i.reshape(-1), j.reshape(-1),k.reshape(-1))
            mapping_from_attribute_to_result.update(result)
            
        results = []
        for attribute in attributes:
            returned_value = mapping_from_attribute_to_result[attribute]
            
            if len(i.shape) == 0:
                value = returned_value[0]
            else:
                value = returned_value.reshape(i.shape)
                
            results.append(value)
            
        return results
        
    def _set_values(self,  indices, attributes, quantities):
        i,j,k = self._to_arrays_of_indices(indices)
        one_dimensional_values = [x.reshape(-1) for x in quantities]
        for setter in self.select_setters_for(attributes):
            setter.set_attribute_values(self, attributes, one_dimensional_values, i.reshape(-1), j.reshape(-1),k.reshape(-1))
     
        
    def _has_key(self, key):
        return key in self.mapping_from_particle_to_index
        
    def _get_keys(self):
        return None 
        
    def __len__(self):
        shape = self.storage_shape()
        return shape[0] * shape[1] * shape[2]
        
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
        return self.writable_attributes


