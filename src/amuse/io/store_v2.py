import h5py
import numpy
import pickle
import os.path

from amuse.units import si
from amuse.units import units
from amuse.units import core
from amuse.units.quantities import is_quantity
from amuse.support import exceptions

from amuse.datamodel import Particles
from amuse.datamodel import Particle
from amuse.datamodel import AttributeStorage
from amuse.datamodel import LinkedArray
from amuse.datamodel import GridPoint
from amuse.datamodel import AbstractSet

from amuse.io import store_v1

import warnings

class HDF5Attribute(object):
    
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
    def new_attribute(cls, name, shape, input, group):
        if is_quantity(input):
            dataset = group.create_dataset(name, shape=shape, dtype=input.number.dtype)
            return HDF5VectorQuantityAttribute(name, dataset, input.unit) 
        elif hasattr(input, 'as_set'):
            subgroup = group.create_group(name)
            group.create_dataset('keys', shape=shape, dtype=input.key.dtype)
            group.create_dataset('masked', shape=shape, dtype=numpy.bool)
            return HDF5LinkedAttribute(name, subgroup) 
        else:
            dtype = numpy.asanyarray(input).dtype
            dataset = group.create_dataset(name, shape=shape, dtype=dtype)
            return HDF5UnitlessAttribute(name, dataset)

    @classmethod
    def load_attribute(cls, name, dataset, loader):
        units_string = dataset.attrs["units"]
        if units_string == "none":
            return HDF5UnitlessAttribute(name, dataset)
        elif units_string == "link":
            return HDF5LinkedAttribute(name, dataset, loader)
        else:
            unit = eval(units_string, core.__dict__)
            return HDF5VectorQuantityAttribute(name, dataset, unit) 

    def get_value(self, index):
        pass

    def remove_indices(self, indices):
        pass

class HDF5VectorQuantityAttribute(HDF5Attribute):
    
    def __init__(self, name, dataset, unit):
        HDF5Attribute.__init__(self, name)
    
        self.dataset = dataset
        self.unit = unit
        
    def get_values(self, indices):
        if indices is None:
            return self.unit.new_quantity(self.dataset[:])
        elif len(self.get_shape()) > 1:
            return self.unit.new_quantity(self.dataset[:][indices])
        else:
            return self.unit.new_quantity(self.dataset[:][indices])
            
    
    def set_values(self, indices, values):
        try:
            self.dataset[indices] = values.value_in(self.unit)
        except AttributeError:
            if not is_quantity(values):
                raise ValueError("Tried to put a non quantity value in a quantity")
            else:
                raise
    
    def get_shape(self):
        return self.dataset.shape
    

    def increase_to_length(self, newlength):
        oldlength = len(self.dataset)
        delta = newlength - len(self.dataset)
        if delta == 0: 
           return
        newshape = list(self.dataset.shape)
        newshape[0] = newlength
        
        values = numpy.empty(shape=self.dataset.shape, dtype=self.dataset.dtype)
        values[:] = self.dataset[:]
    
        parent = self.dataset.parent
        del parent[self.name]
        self.dataset = parent.create_dataset(self.name, shape=newshape, dtype=values.dtype)
        self.dataset[:oldlength] = values[:]

    def get_length(self):
        return len(self.dataset)

    def copy(self):
        #result = type(self)(self.name, self.get_shape(), self.quantity.unit)
        #result.set_values(None, self.get_values(None))
        #return result
        raise NotImplementedError("Copy of an attribute in a hdf5 file is not implemented")

    def get_value(self, index):
        return self.unit.new_quantity(self.dataset[index])

    def remove_indices(self, indices):
        oldlength = len(self.dataset)
        
        values = numpy.empty(shape=self.dataset.shape, dtype=self.dataset.dtype)
        values[:] = self.dataset[:]
        values = numpy.delete(values, indices)
        
        parent = self.dataset.parent
        del parent[self.name]
        self.dataset = parent.create_dataset(self.name, shape=values.shape, dtype=values.dtype)
        self.dataset[:] = values[:]
        
        

    def has_units(self):
        return True

class HDF5UnitlessAttribute(HDF5Attribute):
    
    def __init__(self, name, dataset):
        HDF5Attribute.__init__(self, name)
        
        self.dataset = dataset
        
    def get_values(self, indices):
        if indices is None:
            return self.dataset[:]
        elif len(self.get_shape()) > 1:
            return self.dataset[:][indices]
        else:
            return self.dataset[:][indices]
        
    
    def set_values(self, indices, values):
        self.dataset[indices] = values
    
    def get_shape(self):
        return self.dataset.shape
    
    def increase_to_length(self, newlength):
        oldlength = len(self.dataset)
        delta = newlength - len(self.dataset)
        if delta == 0: 
           return
        newshape = list(self.dataset.shape)
        newshape[0] = newlength
        
        values = numpy.empty(shape=self.dataset.shape, dtype=self.dataset.dtype)
        values[:] = self.dataset[:]
    
        parent = self.dataset.parent
        del parent[self.name]
        self.dataset = parent.create_dataset(self.name, shape=newshape, dtype=values.dtype)
        self.dataset[:oldlength] = values[:]

    def get_length(self):
        return len(self.dataset)

    def copy(self):
        #result = type(self)(self.name, self.get_shape(), self.quantity.unit)
        #result.set_values(None, self.get_values(None))
        #return result
        raise NotImplementedError("Copy of an attribute in hdf5 is not implemented")

    def get_value(self, index):
        return self.dataset[index]

    def remove_indices(self, indices):
        oldlength = len(self.dataset)
        
        values = numpy.empty(shape=self.dataset.shape, dtype=self.dataset.dtype)
        values[:] = self.dataset[:]
        values = numpy.delete(values, indices)
        
        parent = self.dataset.parent
        del parent[self.name]
        self.dataset = parent.create_dataset(self.name, shape=values.shape, dtype=values.dtype)
        self.dataset[:] = values[:]

    def has_units(self):
        return False


class HDF5LinkedAttribute(HDF5Attribute):
    
    def __init__(self, name, group, loader):
        HDF5Attribute.__init__(self, name)
        self.group = group
        self.keys_dataset = group['keys']
        if 'indices' in group:
            self.indices_dataset = group['indices']
        else:
            self.indices_dataset = None
            
        self.kind_dataset = group['kind']
        self.ref_dataset = group['ref']
        self.loader = loader

    def get_values(self, indices):
        kinds = self.kind_dataset[:][indices]
        references = self.ref_dataset[:][indices]
        keys = self.keys_dataset[:][indices]
        if self.indices_dataset:
            grid_indices = self.indices_dataset[:][indices]
        else:
            grid_indices = None
        shape = kinds.shape
        result = LinkedArray(numpy.empty(shape, dtype = numpy.object))
        
        for index, (kind, reference) in enumerate(zip(kinds, references)):
            if   kind == 0:
                result[index] = None
            elif kind == 1:
                
                referenced_group = self.loader.derefence(reference)
                mapping_from_groupid_to_set = self.loader.mapping_from_groupid_to_set
            
                if not referenced_group.id in mapping_from_groupid_to_set:
                    linked_set = self.loader.load_from_group(referenced_group)
                else:
                    linked_set = mapping_from_groupid_to_set[referenced_group.id]
                
                result[index] = linked_set._get_particle(keys[index])
            elif kind == 2:
                referenced_group = self.loader.derefence(reference)
                mapping_from_groupid_to_set = self.loader.mapping_from_groupid_to_set
            
                if not referenced_group.id in mapping_from_groupid_to_set:
                    linked_set = self.loader.load_from_group(referenced_group)
                else:
                    linked_set = mapping_from_groupid_to_set[referenced_group.id]
                    
                result[index] = linked_set._get_gridpoint(tuple(grid_indices[index]))
            elif kind == 3:
                referenced_group = self.loader.derefence(reference)
                mapping_from_groupid_to_set = self.loader.mapping_from_groupid_to_set
            
                if not referenced_group.id in mapping_from_groupid_to_set:
                    linked_set = self.loader.load_from_group(referenced_group)
                else:
                    linked_set = mapping_from_groupid_to_set[referenced_group.id]
                
                result[index] = linked_set
            else:
                raise Exception("unknown link kind")
        return result
        
    
    def set_values(self, indices, values):
        if hasattr(values, 'get_all_keys_in_store'):
            keys = values.get_all_keys_in_store()
            mask = ~values.get_valid_particles_mask()
            self.keys[indices] = keys
            self.masked[indices] = mask
        else:
            if values is None:
                self.keys[indices] = 0
                self.masked[indices] = True
            else:
                for index, key in zip(indices, values):
                    if key is None:
                        self.keys[indices] = 0
                        self.masked[indices] = True
                    else:
                        self.keys[index] = key
                        self.masked[indices] = False
    
    def get_length(self):
        return len(self.keys)
                
    
    def increase_to_length(self, newlength):
        oldlength = len(self.keys)
        delta = newlength - len(self.keys)
        if delta == 0: 
           return
        newshape = list(self.keys.shape)
        newshape[0] = newlength
        values = numpy.empty(shape=self.keys.shape, dtype=self.keys.dtype)
        values[:] = self.keys[:]
    
        parent = self.group
        del parent['keys']
        self.keys = parent.create_dataset('keys', shape=newshape, dtype=values.dtype)
        self.keys[:oldlength] = values[:]
        
        values = numpy.empty(shape=self.masked.shape, dtype=self.masked.dtype)
        values[:] = self.masked[:]
        parent = self.group
        del parent['masked']
        self.masked = parent.create_dataset('masked', shape=newshape, dtype=values.dtype)
        self.masked[:oldlength] = values[:]
        

    def get_shape(self):
        return self.keys.shape

    def copy(self):
        #result = type(self)(self.name, self.get_shape(), self.linked_set)
        #result.set_values(None, self.get_values(None))
        #return result
        raise NotImplementedError("Copy of an attribute in hdf5 is not implemented")

    def get_value(self, index):
        #key = self.values[index]
        #if key is ma.masked:
        #    return None
        #else:
        #    return self.linked_set._get_particle_unsave(key)
        raise NotImplementedError("Copy of an attribute in hdf5 is not implemented")

    def remove_indices(self, indices):
        if 1:
            
            raise NotImplementedError("Copy of an attribute in hdf5 is not implemented")
        oldlength = len(self.dataset)
        
        values = numpy.empty(shape=self.dataset.shape, dtype=self.dataset.dtype)
        values[:] = self.dataset[:]
        values = numpy.delete(values, indices)
        
        parent = self.dataset.parent
        del parent[self.name]
        self.dataset = parent.create_dataset(self.name, shape=values.shape, dtype=values.dtype)
        self.dataset[:] = values[:]

    def has_units(self):
        return False


class HDF5AttributeStorage(AttributeStorage):

    def __init__(self, keys, hdfgroup, loader):
        self.hdfgroup = hdfgroup
        self.attributesgroup = self.hdfgroup["attributes"]
        self.number_of_particles = self.hdfgroup.attrs["number_of_particles"]
        self.particle_keys = keys
        self.loader = loader
        self.mapping_from_particle_to_index = self.new_index()
        
    def can_extend_attributes(self):
        return True
        
    def __len__(self):
        return len(self.particle_keys)
        
    def new_index(self):
        result = {}
        index = 0
        for particle_key in self.particle_keys:
            result[particle_key] = index
            index += 1
        return result
        
    def get_indices_of(self, keys):
        if keys is None:
            return numpy.arange(0,len(self.particle_keys))
            
        mapping_from_particle_to_index = self.mapping_from_particle_to_index 
        result = numpy.zeros(len(keys),dtype='int32')
        
        index = 0
        for index, particle_key in enumerate(keys):
            result[index] = mapping_from_particle_to_index[particle_key]
            index += 1
        return result
        
    def get_defined_attribute_names(self):
        return self.attributesgroup.keys()
    
    def get_values_in_store(self, indices, attributes):
        results = []
        for attribute in attributes:
            dataset = HDF5Attribute.load_attribute(
                attribute,
                self.attributesgroup[attribute],
                self.loader
            )
            selected_values = dataset.get_values(indices)
            results.append(selected_values)
        
        return results
        
    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_to_index
    
    def get_all_keys_in_store(self):
        return self.particle_keys
        
    def set_values_in_store(self, indices, attributes, quantities):
            
        for attribute, quantity in zip(attributes, quantities):
            if attribute in self.attributesgroup:
                dataset = HDF5Attribute.load_attribute(
                    attribute,
                    self.attributesgroup[attribute],
                    self.loader
                )
            else:
                dataset = HDF5Attribute.new_attribute(
                    attribute,
                    len(self.particle_keys),
                    quantity,
                    self.attributesgroup
                )
                
            bools = numpy.zeros(dataset.get_shape(), dtype='bool')
            bools[indices] = True
            dataset.set_values(bools, quantity)

    def get_all_indices_in_store(self):
        return numpy.arange(len(self.particle_keys))

class HDF5GridAttributeStorage(AttributeStorage):

    def __init__(self, shape, hdfgroup, loader):
        self.hdfgroup = hdfgroup
        self.attributesgroup = self.hdfgroup["attributes"]
        self.shape = shape
        self.loader = loader
        
    def storage_shape(self):
        return self.shape
    
    def add_particles_to_store(self, keys, attributes = [], quantities = []):
        raise exceptions.AmuseException("adding points to the grid is not implemented")
            
    def remove_particles_from_store(self, keys):
        raise exceptions.AmuseException("removing points from the grid is not implemented")
        
    def __len__(self):
        return numpy.prod(self.shape)
        
    def get_unit_of(self, attribute):
        dataset = self.attributesgroup[attribute]
        return eval(dataset.attrs["units"], core.__dict__) 
        
    def get_defined_attribute_names(self):
        return self.attributesgroup.keys()
        
    def get_values_in_store(self, indices, attributes):
            
        results = []
        for attribute in attributes:
            dataset = HDF5Attribute.load_attribute(
                attribute,
                self.attributesgroup[attribute],
                self.loader
            )
            selected_values = dataset.get_values(indices)
            results.append(selected_values)
        
        return results
        
    def has_key_in_store(self, key):
        return False
    
    def get_all_keys_in_store(self):
        return None
        
    def set_values_in_store(self, indices, attributes, quantities):
        for attribute, quantity in zip(attributes, quantities):
            if attribute in self.attributesgroup:
                dataset = HDF5Attribute.load_attribute(
                    attribute,
                    self.attributesgroup[attribute],
                    self.loader
                )
            else:
                dataset = HDF5Attribute.new_attribute(
                    attribute,
                    self.shape,
                    quantity,
                    self.attributesgroup
                )
                
            bools = numpy.zeros(dataset.get_shape(), dtype='bool')
            bools[indices] = True
            dataset.set_values(bools, quantity)
        
        
class StoreHDF(object):
    INFO_GROUP_NAME = 'AMUSE_INF'
    DATA_GROUP_NAME = 'data'
    
    def __init__(self, filename, append_to_file=True, open_for_writing = True, copy_history = False, return_working_copy = False):
        if not append_to_file and open_for_writing and os.path.exists(filename):
            os.remove(filename)
            
        if append_to_file:
            if open_for_writing:
                self.hdf5file = h5py.File(filename,'a')
            else:
                if os.access(filename, os.W_OK):
                    self.hdf5file = h5py.File(filename,'a')
                else:
                    self.hdf5file = h5py.File(filename,'r')
        else:
            if open_for_writing:
                self.hdf5file = h5py.File(filename,'w')
            else:
                self.hdf5file = h5py.File(filename,'r')
        
        self.copy_history = copy_history
        self.return_working_copy = return_working_copy
        self.mapping_from_groupid_to_set = {}
        
        warnings.warn("amuse hdf storage version 2.0 is still in development, do not use it for production scripts")
        
    def is_correct_version(self):
        if len(self.hdf5file) == 0 or self.INFO_GROUP_NAME in self.hdf5file:
            return True
        else:
            return False
    
    def get_version(self):
        return "2.0"
        
    def store_sets(self, sets, names, extra_attributes = {}):
        info_group = self.info_group()
        info_group.attrs["version"] = self.get_version()
    
        mapping_from_setid_to_group = {}
        links_to_resolve = []
        
        for x, name in zip(sets, names):
            if hasattr(x, 'shape') and len(x.shape) > 1:
                self.store_grid(x, extra_attributes, parent = self.named_group(name))
            else:   
                self.store_particles(x, extra_attributes, self.named_group(name), mapping_from_setid_to_group, links_to_resolve)
        
        while len(links_to_resolve) > 0:
            sets_to_store, links_to_resolve = self.resolve_links(
                mapping_from_setid_to_group,
                links_to_resolve
            )
            for group_to_store_under, x in sets_to_store:
                if hasattr(x, 'shape') and len(x.shape) > 1:
                    self.store_grid(x, {}, group_to_store_under, mapping_from_setid_to_group, links_to_resolve)
                else:   
                    self.store_particles(x, {}, group_to_store_under, mapping_from_setid_to_group, links_to_resolve)
        
        
    def store(self, container, extra_attributes = {}):
        info_group = self.info_group()
        info_group.attrs["version"] = self.get_version()
        
        mapping_from_setid_to_group = {}
        links_to_resolve = []
            
        if hasattr(container, 'keys') and not hasattr(container, 'as_set'):
            self.store_sets(
                container.values(),
                container.keys(),
                extra_attributes
            )
        if hasattr(container, 'shape') and len(container.shape) > 1:
            self.store_grid(
                container, 
                extra_attributes,
                None,
            )
        else:      
            self.store_particles(
                container, 
                extra_attributes, 
                None, 
                mapping_from_setid_to_group, 
                links_to_resolve
            )
        
        while len(links_to_resolve) > 0:
            sets_to_store, links_to_resolve = self.resolve_links(
                mapping_from_setid_to_group,
                links_to_resolve
            )
            
            for group_to_store_under, x in sets_to_store:
                if hasattr(x, 'shape') and len(x.shape) > 1:
                    self.store_grid(x, {}, group_to_store_under, mapping_from_setid_to_group, links_to_resolve)
                else:   
                    self.store_particles(x, {}, group_to_store_under, mapping_from_setid_to_group, links_to_resolve)
        
    def store_particles(self, particles, extra_attributes = {}, parent=None, mapping_from_setid_to_group = {}, links = []):
        if parent is None:
            parent = self.data_group()
            
        group = self.new_version(parent)
        group.attrs["type"] = 'particles'
        self.mapping_from_groupid_to_set[group.id] = particles._original_set()
        
        
        group.attrs["number_of_particles"] = len(particles)
        group.attrs["class_of_the_particles"] = pickle.dumps(particles._factory_for_new_collection())
            
        keys = particles.get_all_keys_in_store()
        dataset = group.create_dataset("keys", data=keys)
        self.hdf5file.flush()
        self.store_collection_attributes(particles, group, extra_attributes)
        self.store_values(particles, group, links)
            
        mapping_from_setid_to_group[id(particles._original_set())] = group
        
        self.hdf5file.flush()
        
    
    def store_grid(self, grid, extra_attributes = {}, parent=None, mapping_from_setid_to_group = {}, links = []):
        if parent is None:
            parent = self.data_group()
            
        group = self.new_version(parent)
        
        group.attrs["type"] = 'grid'
        group.attrs["class_of_the_container"] = pickle.dumps(grid._factory_for_new_collection())
        group.create_dataset("shape", data=numpy.asarray(grid.shape))
    
        self.store_collection_attributes(grid, group, extra_attributes)
        self.store_values(grid, group)
        
        mapping_from_setid_to_group[id(grid._original_set())] = group
        
        self.hdf5file.flush()
        
    
    def resolve_links(self, mapping_from_setid_to_group, links):
        sets_to_store = []
        seen_sets_by_id = set([])
        links_to_resolve = []
        for group, index, ref_dataset, linked_set in links:
            if id(linked_set) in mapping_from_setid_to_group:
                stored_group = mapping_from_setid_to_group[id(linked_set)]
                ref_dataset[index] = stored_group.ref
            else:
                if not id(linked_set) in seen_sets_by_id:
                    sets_to_store.append([group, linked_set])
                    seen_sets_by_id.add(id(linked_set))
                links_to_resolve.append([group, index, ref_dataset, linked_set])
        return sets_to_store, links_to_resolve
        
    def store_values(self, container, group, links = []):
        attributes_group = group.create_group("attributes")
        
        all_values = container.get_values_in_store(None, container.get_attribute_names_defined_in_store())
        for attribute, quantity in zip(container.get_attribute_names_defined_in_store(), all_values):
            if is_quantity(quantity):
                value = quantity.value_in(quantity.unit)
                dataset = attributes_group.create_dataset(attribute, data=value)
                dataset.attrs["units"] = quantity.unit.to_simple_form().reference_string()
            elif isinstance(quantity, LinkedArray):
                subgroup = attributes_group.create_group(attribute)
                shape = quantity.shape
                kind_array = numpy.zeros(shape, dtype = numpy.int16)
                ref_dtype = h5py.special_dtype(ref=h5py.Reference)
                ref_array = numpy.empty(shape, dtype = ref_dtype)
                ref_dataset = subgroup.create_dataset('ref', data=ref_array)
                key_array = numpy.zeros(shape, dtype = numpy.uint64)
                
                max_len_grid_indices = 0
                for index, object in enumerate(quantity):
                    if isinstance(object, GridPoint):
                        max_len_grid_indices = max(len(object.index), max_len_grid_indices)
                
                if max_len_grid_indices > 0:
                    indices_shape = list(shape)
                    indices_shape.append(max_len_grid_indices)
                    indices_array = numpy.zeros(indices_shape, dtype = numpy.uint64)
                    
                for index, object in enumerate(quantity):
                    if object is None:
                        kind_array[index] = 0
                        key_array[index] = 0
                    elif isinstance(object, Particle):
                        kind_array[index] = 1
                        key_array[index] = object.key
                        links.append([group, index, ref_dataset, object.get_containing_set()])
                    elif isinstance(object, GridPoint):
                        kind_array[index] = 2
                        key_array[index] = 0
                        grid_index = object.index
                        if len(grid_index) < max_len_grid_indices:
                            grid_index = list(grid_index)
                            for _ in range(max_len_grid_indices-len(grid_index)):
                                grid_index.append(0)
                            
                        indices_array[index] = grid_index
                        links.append([group, index, ref_dataset, object.get_containing_set()])
                    elif isinstance(object, AbstractSet):
                        kind_array[index] = 3
                        key_array[index] = 0
                        links.append([group, index,  ref_dataset, object._original_set()])
                    else:
                        raise Exception("unsupported type")
                if max_len_grid_indices > 0:
                    subgroup.create_dataset('indices', data=indices_array)
                    
                subgroup.create_dataset('keys', data=key_array)
                subgroup.create_dataset('kind', data=kind_array)
                subgroup.attrs["units"] = "link"
            else:
                dataset = attributes_group.create_dataset(attribute, data=quantity)
                dataset.attrs["units"] = "none"
                
            
    
    def store_timestamp(self, container, group):
        quantity = container.get_timestamp()
        if not quantity is None:
            group.attrs["timestamp"] = quantity.value_in(quantity.unit)
            group.attrs["timestamp_unit"] = quantity.unit.reference_string()
    
    
    
    def store_collection_attributes(self, container, group, extra_attributes):
        collection_attributes = container.collection_attributes.__getstate__()
        
        arguments_and_attributes = {}
        arguments_and_attributes.update(collection_attributes)
        arguments_and_attributes.update(extra_attributes)
        
        for name, quantity in arguments_and_attributes.iteritems():
            if quantity is None:
                continue 
            if is_quantity(quantity):
                group.attrs[name] = quantity.value_in(quantity.unit)
                group.attrs[name+"_unit"] = quantity.unit.reference_string()
            else:
                group.attrs[name] = quantity
                group.attrs[name+"_unit"] = "none"
            
    def load_collection_attributes(self, container, group):
        names = group.attrs.keys()
        attributenames = [x for x in names if x + '_unit' in group.attrs]
        for name in attributenames:
            unit_string = group.attrs[name+"_unit"]
            if unit_string == 'none':
                quantity = group.attrs[name]
            else:
                unit = eval(group.attrs[name+"_unit"], core.__dict__) 
                quantity = unit.new_quantity(group.attrs[name])
            setattr(container.collection_attributes, name, quantity)
                
    def load(self):
        if len(self.data_group()) > 0:
            return self.load_container(self.data_group())
        else:
            pass
                
    def load_sets(self, names):
        result = []
        for x in names:
            set = self.load_container(self.named_group(x))
            result.append(set)
        return result
            
    def load_particles_from_group(self, group):
        try:
            class_of_the_container = pickle.loads(group.attrs["class_of_the_particles"])
        except:
            class_of_the_container = Particles
            
        dataset = group["keys"]
        keys = numpy.ndarray(len(dataset), dtype = dataset.dtype)
        dataset.read_direct(keys)
        
        particles = class_of_the_container(is_working_copy = False)
        particles._private.attribute_storage = HDF5AttributeStorage(keys, group, self)
        self.load_collection_attributes(particles, group)
        
        self.mapping_from_groupid_to_set[group.id] = particles
        
        return particles
        
        
    def load_grid_from_group(self, group):
        try:
            class_of_the_container = pickle.loads(group.attrs["class_of_the_container"])
        except:
            class_of_the_container = Grids
            
        shape = tuple(group["shape"])
        
        container = class_of_the_container()
        container._private.attribute_storage = HDF5GridAttributeStorage(shape, group, self)
        self.load_collection_attributes(container, group)
        
        self.mapping_from_groupid_to_set[group.id] = container
        return container
    
    def load_from_group(self, group):
        container_type = group.attrs['type']
        
        if container_type == 'particles':
            return self.load_particles_from_group(group)
        elif container_type == 'grid':
            return self.load_grid_from_group(group)
        else:
            raise Exception('unknown container type in file {0}'.format(container_type))
        
        
    def load_particles(self, container_group = None):
        if container_group is None:
            container_group = self.data_group()
            
        return self.load_container(container_group)
        
    
    def load_grid(self, container_group = None):
        if container_group is None:
            container_group = self.data_group()
            
        return self.load_container(container_group)
        
    def load_container(self, container_group):
        number_of_saved_containers= len(container_group)
        all_containers = [None] * number_of_saved_containers
        for group_index in container_group.keys():
            group = container_group[group_index]
            container = self.load_from_group(group)
            if self.copy_history:
                container = container.copy()
            all_containers[int(group_index) - 1] = container
            
        previous = None
        for x in all_containers:
            x._private.previous = previous
            previous = x
            
        last = all_containers[-1]
        return last
    
    def derefence(self, reference):
        return self.hdf5file[reference]
        
    def new_version(self, master_group):
        index = len(master_group)
        name = format(index + 1,"010d")
        return master_group.create_group(name)
        
    def info_group(self):
        return self.named_group(self.INFO_GROUP_NAME)
        
    def data_group(self):
        return self.named_group(self.DATA_GROUP_NAME)
        
    def named_group(self, name):
        return self.hdf5file.require_group(name)
        
    def close(self):
        self.hdf5file.flush()
        self.hdf5file.close()
