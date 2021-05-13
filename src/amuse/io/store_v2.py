try:
    import h5py
except ImportError as ex:
    import warnings
    warnings.warn("could not load h5py, hdf5 files not supported", ImportWarning)
    h5py = None
    
import numpy
import pickle
import os.path
import sys

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

import logging
logger = logging.getLogger(__name__)

def pickle_to_string(value):
    return numpy.void(pickle.dumps(value, protocol=0))
        
def unpickle_from_string(value):
    return pickle.loads(value, encoding='bytes')
    
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
            if not hasattr(shape, '__iter__'): 
                shape = shape,
            dtype = numpy.asanyarray(input.number).dtype
            dataset = group.create_dataset(name, shape=shape, dtype=dtype)
            dataset.attrs["units"] = input.unit.to_simple_form().reference_string().encode('ascii')
            return HDF5VectorQuantityAttribute(name, dataset, input.unit)                                     
        elif hasattr(input, 'as_set'):
            raise Exception("adding a linked attribute to a set stored in a HDF file is not supported, alternative is to copy the set and save it")
            subgroup = group.create_group(name)
            group.create_dataset('keys', shape=shape, dtype=input.key.dtype)
            group.create_dataset('masked', shape=shape, dtype=numpy.bool)
            return HDF5LinkedAttribute(name, subgroup)                                     
        else:
            dtype = numpy.asanyarray(input).dtype
            if dtype.kind == 'U':
                new_dtype = numpy.dtype('S' + dtype.itemsize * 4)
                dataset = group.create_dataset(name, shape=shape, dtype=dtype)
                dataset.attrs["units"] = "UNICODE".encode('ascii')
                return HDF5UnicodeAttribute(name, dataset)
            else:
                if not hasattr(shape, '__iter__'): 
                    shape = shape,
                dataset = group.create_dataset(name, shape=shape, dtype=dtype)
                dataset.attrs["units"] = "none".encode('ascii')
                return HDF5UnitlessAttribute(name, dataset)



    @classmethod
    def load_attribute(cls, name, dataset, loader):
        units_string = dataset.attrs["units"] if isinstance(dataset.attrs["units"], str) else dataset.attrs["units"].decode("ascii") 
        if units_string == "none":
            return HDF5UnitlessAttribute(name, dataset)
        elif units_string == "link":
            return HDF5LinkedAttribute(name, dataset, loader)
        elif units_string == "UNICODE":
            return HDF5UnicodeAttribute(name, dataset)
        else:
            try:
                unit = eval(units_string, core.__dict__)
            except:
                # These are some strange things found in the
                # unit string field of a couple of data files
                # I'm unsure if this is a data corruption
                # or some problem in printing units...
                # same error was found in multiple files
                # pointing to a real error
                # this code is for existing files
                if '*0(' in units_string:
                    updated = units_string.replace('*0(','* (')
                    unit = eval(updated, core.__dict__)
                elif 'vystem.' in units_string:
                    updated = units_string.replace('vystem.','system.')
                    unit = eval(updated, core.__dict__)
                else:
                    print(units_string)
                    raise

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
        if indices is None: 
            indices=slice(None)
        kinds = self.kind_dataset[:][indices]
        references = self.ref_dataset[:][indices]
        keys = self.keys_dataset[:][indices]
        shape = kinds.shape
        if self.indices_dataset:
            grid_indices = self.indices_dataset[:][indices]
        else:
            grid_indices = numpy.zeros(shape)
        result = LinkedArray(numpy.empty(shape, dtype = numpy.object))
        if len(shape) == 0:
            # we have one unique value, happens with grids
            result = self.convert_to_object(kinds, references, keys, grid_indices)
        else:
            for index in numpy.ndindex(*shape):
                if len(index) == 1:
                    index = index[0]
                reference = references[index]
                kind = kinds[index]
                result[index] = self.convert_to_object(kind,reference, keys[index], grid_indices[index])
            
        return result
        
    
    def convert_to_object(self, kind, reference, key, grid_index):
        if   kind == 0:
            return None
        elif kind == 1:
            referenced_group = self.loader.derefence(reference)
            
            mapping_from_groupid_to_set = self.loader.mapping_from_groupid_to_set
        
            if not referenced_group.id in mapping_from_groupid_to_set:
                linked_set = self.loader.load_from_group(referenced_group)
            else:
                linked_set = mapping_from_groupid_to_set[referenced_group.id]
            
            return linked_set._get_particle(key)
        elif kind == 2:
            referenced_group = self.loader.derefence(reference)
            mapping_from_groupid_to_set = self.loader.mapping_from_groupid_to_set
        
            if not referenced_group.id in mapping_from_groupid_to_set:
                linked_set = self.loader.load_from_group(referenced_group)
            else:
                linked_set = mapping_from_groupid_to_set[referenced_group.id]
                
            return linked_set._get_gridpoint(tuple(grid_index))
        elif kind == 3:
            referenced_group = self.loader.derefence(reference)
            mapping_from_groupid_to_set = self.loader.mapping_from_groupid_to_set
        
            if not referenced_group.id in mapping_from_groupid_to_set:
                linked_set = self.loader.load_from_group(referenced_group)
            else:
                linked_set = mapping_from_groupid_to_set[referenced_group.id]
            
            return linked_set
        else:
            raise Exception("unknown link kind")
            
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
        return list(self.attributesgroup.keys())
        
    def get_defined_settable_attribute_names(self):
        return self.get_defined_attribute_names()
    
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
        decoded = dataset.attrs["units"] if isinstance(dataset.attrs["units"], str) else dataset.attrs["units"].decode("ascii")
        return eval(decoded, core.__dict__) 
        
    def get_defined_attribute_names(self):
        return list(self.attributesgroup.keys())
        
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
        return Ellipsis
        
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
    
class UneresolvedItemInArrayLink(object):

    def __init__(self, group, index, dataset_to_resolve, linked_set):
        self.group = group
        self.index = index
        self.dataset_to_resolve = dataset_to_resolve
        self.linked_set = linked_set
        self.resolved = False
        
    def resolve(self, mapping_from_setid_to_group):
        if id(self.linked_set) in mapping_from_setid_to_group:
            stored_group = mapping_from_setid_to_group[id(self.linked_set)]
            self.dataset_to_resolve[self.index] = stored_group.ref
            self.resolved = True
        else:
            self.resolved = False
    
    def is_resolved(self):
        return self.resolved
        
        
class UneresolvedAttributeLink(object):

    def __init__(self, group, name, linked_set):
        self.group = group
        self.name = name
        self.linked_set = linked_set
        self.resolved = False
        
    def resolve(self, mapping_from_setid_to_group):
        if id(self.linked_set) in mapping_from_setid_to_group:
            stored_group = mapping_from_setid_to_group[id(self.linked_set)]
            self.group.attrs[self.name] = stored_group.ref
            self.resolved = True
        else:
            self.resolved = False
    
    def is_resolved(self):
        return self.resolved
        
class StoreHDF(object):
    INFO_GROUP_NAME = 'AMUSE_INF'
    DATA_GROUP_NAME = 'data'
    
    def __init__(self, filename, append_to_file=True, open_for_writing = True, copy_history = False, overwrite_file=False):
        if h5py is None:
            raise exceptions.AmuseException("h5py module not available, cannot use hdf5 files")
        
        logger.info("opening {0} with options {1} {2} {3} {4}".format(filename, append_to_file, open_for_writing, copy_history,overwrite_file))
            
        if not append_to_file and open_for_writing:
            if os.path.exists(filename):
                if overwrite_file:
                    os.remove(filename)
                else:
                    raise FileExistsError("Opening file for write with overwrite_file is False but file {0} exists".format(filename))

        if append_to_file:
            if os.access(filename, os.F_OK) and not os.access(filename, os.W_OK):
                   raise Exception("Opening file for append but file {0} is not writeable".format(filename))
            self.hdf5file = h5py.File(filename,'a')
        elif open_for_writing:
            self.hdf5file = h5py.File(filename,'w')
        else:
            self.hdf5file = h5py.File(filename,'r')
        
        self.copy_history = copy_history
        self.mapping_from_groupid_to_set = {}
        
        #warnings.warn("amuse hdf storage version 2.0 is still in development, do not use it for production scripts")
        
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
            if hasattr(x, 'shape'):
                self.store_grid(x, extra_attributes, parent = self.named_group(name))
            else:   
                self.store_particles(x, extra_attributes, self.named_group(name), mapping_from_setid_to_group, links_to_resolve)
        
        while len(links_to_resolve) > 0:
            sets_to_store, links_to_resolve = self.resolve_links(
                mapping_from_setid_to_group,
                links_to_resolve
            )
            for group_to_store_under, x in sets_to_store:
                if hasattr(x, 'shape'):
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
                list(container.values()),
                list(container.keys()),
                extra_attributes
            )
        if hasattr(container, 'shape'):
            self.store_grid(
                container, 
                extra_attributes,
                None,
                mapping_from_setid_to_group, 
                links_to_resolve
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
                if hasattr(x, 'shape'):
                    self.store_grid(x, {}, group_to_store_under, mapping_from_setid_to_group, links_to_resolve)
                else:   
                    self.store_particles(x, {}, group_to_store_under, mapping_from_setid_to_group, links_to_resolve)
        
    def store_particles(self, particles, extra_attributes = {}, parent=None, mapping_from_setid_to_group = {}, links = []):
        if parent is None:
            parent = self.data_group()
            
        group = self.new_version(parent)
        group.attrs["type"] = 'particles'.encode("ascii")
        self.mapping_from_groupid_to_set[group.id] = particles
        
        
        group.attrs["number_of_particles"] = len(particles)
        group.attrs["class_of_the_particles"] = pickle_to_string(particles._factory_for_new_collection())
            
        keys = particles.get_all_keys_in_store()
        dataset = group.create_dataset("keys", data=keys)
        self.hdf5file.flush()
        self.store_collection_attributes(particles, group, extra_attributes, links)
        self.store_values(particles, group, links)
            
        mapping_from_setid_to_group[id(particles)] = group
        
        self.hdf5file.flush()
        
    
    def store_grid(self, grid, extra_attributes = {}, parent=None, mapping_from_setid_to_group = {}, links = []):
        if parent is None:
            parent = self.data_group()
            
        group = self.new_version(parent)
        
        group.attrs["type"] = 'grid'.encode("ascii")
        group.attrs["class_of_the_container"] = pickle_to_string(grid._factory_for_new_collection())
        group.create_dataset("shape", data=numpy.asarray(grid.shape))
    
        self.store_collection_attributes(grid, group, extra_attributes, links)
        self.store_values(grid, group, links)
        
        mapping_from_setid_to_group[id(grid)] = group
        
        self.hdf5file.flush()
        
    
    def resolve_links(self, mapping_from_setid_to_group, links):
        sets_to_store = []
        seen_sets_by_id = set([])
        links_to_resolve = []
        for unresolved_link in links:
            unresolved_link.resolve(mapping_from_setid_to_group)
            if not unresolved_link.is_resolved():
                if not id(unresolved_link.linked_set) in seen_sets_by_id:
                    sets_to_store.append([unresolved_link.group, unresolved_link.linked_set])
                    seen_sets_by_id.add(id(unresolved_link.linked_set))
                links_to_resolve.append(unresolved_link)
        return sets_to_store, links_to_resolve
        
    def store_values(self, container, group, links = []):
        attributes_group = group.create_group("attributes")
        
        all_values = container.get_values_in_store(Ellipsis, container.get_attribute_names_defined_in_store())
        for attribute, quantity in zip(container.get_attribute_names_defined_in_store(), all_values):
            if is_quantity(quantity):
                value = quantity.value_in(quantity.unit)
                dataset = attributes_group.create_dataset(attribute, data=value)
                dataset.attrs["units"] = quantity.unit.to_simple_form().reference_string().encode("ascii")
            elif isinstance(quantity, LinkedArray):
                self.store_linked_array(attribute, attributes_group, quantity, group, links)
            else:
                dtype = numpy.asanyarray(quantity).dtype
                if dtype.kind == 'U':
                    dataset = attributes_group.create_dataset(attribute, data=numpy.char.encode(quantity,  'UTF-32BE'))
                    dataset.attrs["units"] = "UNICODE".encode('ascii')
                else:
                    dataset = attributes_group.create_dataset(attribute, data=quantity)
                    dataset.attrs["units"] = "none".encode('ascii')
                
    

    def store_linked_array(self, attribute, attributes_group, quantity, group, links):
        subgroup = attributes_group.create_group(attribute)
        shape = quantity.shape
        kind_array = numpy.zeros(shape, dtype = numpy.int16)
        ref_dtype = h5py.special_dtype(ref=h5py.Reference)
        ref_array = numpy.empty(shape, dtype = ref_dtype)
        ref_dataset = subgroup.create_dataset('ref', data=ref_array)
        key_array = numpy.zeros(shape, dtype = numpy.uint64)
        
        max_len_grid_indices = 0
        for index in numpy.ndindex(*shape):
            object = quantity[index]
            if isinstance(object, GridPoint):
                max_len_grid_indices = max(len(object.index), max_len_grid_indices)
        
        if max_len_grid_indices > 0:
            indices_shape = list(shape)
            indices_shape.append(max_len_grid_indices)
            indices_array = numpy.zeros(indices_shape, dtype = numpy.uint64)
            
        for index in numpy.ndindex(*shape):
            object = quantity[index]
            if object is None:
                kind_array[index] = 0
                key_array[index] = 0
            elif isinstance(object, Particle):
                kind_array[index] = 1
                key_array[index] = object.key
                links.append(UneresolvedItemInArrayLink(
                    group, 
                    index, 
                    ref_dataset, 
                    object.get_containing_set()
                ))
            elif isinstance(object, GridPoint):
                kind_array[index] = 2
                key_array[index] = 0
                grid_index = object.index
                if len(grid_index) < max_len_grid_indices:
                    grid_index = list(grid_index)
                    for _ in range(max_len_grid_indices-len(grid_index)):
                        grid_index.append(0)
                    
                indices_array[index] = grid_index
                links.append(UneresolvedItemInArrayLink(
                    group, 
                    index, 
                    ref_dataset, 
                    object.get_containing_set()
                ))
            elif isinstance(object, AbstractSet):
                kind_array[index] = 3
                key_array[index] = 0
                links.append(UneresolvedItemInArrayLink(
                    group, 
                    index, 
                    ref_dataset, 
                    object#.copy() #._original_set()
                ))
            else:
                raise Exception("unsupported type {0}".format(type(object)))
        if max_len_grid_indices > 0:
            subgroup.create_dataset('indices', data=indices_array)
            
        subgroup.create_dataset('keys', data=key_array)
        subgroup.create_dataset('kind', data=kind_array)
        subgroup.attrs["units"] = "link".encode("ascii")
    

    def store_timestamp(self, container, group):
        quantity = container.get_timestamp()
        if not quantity is None:
            group.attrs["timestamp"] = quantity.value_in(quantity.unit)
            group.attrs["timestamp_unit"] = quantity.unit.reference_string().encode("ascii")
    
    
    
    def store_collection_attributes(self, container, group, extra_attributes, links):
        collection_attributes = container.collection_attributes.__getstate__()
        arguments_and_attributes = {}
        arguments_and_attributes.update(collection_attributes)
        arguments_and_attributes.update(extra_attributes)
        ref_dtype = h5py.special_dtype(ref=h5py.Reference)
        for name, quantity in arguments_and_attributes.items():
            if quantity is None:
                continue 
            if is_quantity(quantity):
                group.attrs[name] = quantity.value_in(quantity.unit)
                group.attrs[name+"_unit"] = quantity.unit.reference_string().encode("ascii")
            elif isinstance(quantity, Particle):
                #group.attrs[name] = ref_dtype(None)
                group.attrs[name+"_key"] = quantity.key
                group.attrs[name+"_unit"] = "particle".encode("ascii")
                links.append(UneresolvedAttributeLink(
                    group,
                    name,
                    quantity.get_containing_set()
                ))
            elif isinstance(quantity, GridPoint):
                #group.attrs[name] = ref_dtype(None)
                group.attrs[name+"_index"] = quantity.index
                group.attrs[name+"_unit"] = "gridpoint".encode("ascii")
                links.append(UneresolvedAttributeLink(
                    group, 
                    name, 
                    quantity.get_containing_set()
                ))
            elif isinstance(quantity, AbstractSet):
                #group.attrs[name] = ref_dtype(None)
                group.attrs[name+"_unit"] = "set".encode("ascii")
                links.append(UneresolvedAttributeLink(
                    group, 
                    name,
                    quantity._original_set()
                ))
            else:
                group.attrs[name] = quantity
                group.attrs[name+"_unit"] = "none".encode("ascii")
            
    def load_collection_attributes(self, container, group):
        names = group.attrs.keys()
        attributenames = [x for x in names if x + '_unit' in group.attrs]
        for name in attributenames:
            unit_string = group.attrs[name+"_unit"] if isinstance(group.attrs[name+"_unit"],str) else group.attrs[name+"_unit"].decode("ascii")
            if unit_string == 'none':
                quantity = group.attrs[name]
            elif unit_string == 'particle':
                reference = group.attrs[name]
                set = self.get_set_from_reference(reference)
                key = group.attrs[name+'_key']
                quantity = set._get_particle_unsave(key)
            elif unit_string == 'gridpoint':
                reference = group.attrs[name]
                set = self.get_set_from_reference(reference)
                index = group.attrs[name+'_index']
                quantity = set._get_gridpoint(tuple(index))
            elif unit_string == 'set':
                reference = group.attrs[name]
                quantity = self.get_set_from_reference(reference)
            else:
                unit = eval(unit_string, core.__dict__) 
                quantity = unit.new_quantity(group.attrs[name])
            setattr(container.collection_attributes, name, quantity)
                
    def load(self):
        if not self.data_group(False) is None and len(self.data_group(False)) > 0:
            return self.load_container(self.data_group())
        else:
            result = {}
            for x in self.hdf5file.keys():
                if x == 'AMUSE_INF':
                    continue
                result[x] = self.load_container(self.named_group(x))
            return result
                
    def load_sets(self, names):
        result = []
        for x in names:
            set = self.load_container(self.named_group(x))
            result.append(set)
        return result
            
    def load_particles_from_group(self, group):
        try:
            class_of_the_container = unpickle_from_string(group.attrs["class_of_the_particles"])
        except:
            class_of_the_container = Particles
            
        dataset = group["keys"]
        keys = numpy.ndarray(len(dataset), dtype = dataset.dtype)
        if len(keys) == 0:
            particles = class_of_the_container(is_working_copy = False)
            self.mapping_from_groupid_to_set[group.id] = particles
        else:
            dataset.read_direct(keys)
        
            particles = class_of_the_container(is_working_copy = False)
            particles._private.attribute_storage = HDF5AttributeStorage(keys, group, self)
       
            self.mapping_from_groupid_to_set[group.id] = particles
            self.load_collection_attributes(particles, group)
        
        
        return particles
        
    def load_grid_from_group(self, group):
        try:
            class_of_the_container = unpickle_from_string(group.attrs["class_of_the_container"])
        except:
            class_of_the_container = Grids
            
        shape = tuple(group["shape"])
        
        container = class_of_the_container()
        container._private.attribute_storage = HDF5GridAttributeStorage(shape, group, self)
        self.mapping_from_groupid_to_set[group.id] = container
        self.load_collection_attributes(container, group)
        
        return container
    
    def load_from_group(self, group):
        container_type = group.attrs['type'] if isinstance(group.attrs['type'], str) else group.attrs['type'].decode('ascii') 
        
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
        
        if self.copy_history:
            copy_of_last = last.copy()
            copy_of_last._private.previous = last
            return copy_of_last
        else:
            return last
    
    def get_set_from_reference(self, reference):
        referenced_group = self.derefence(reference)
        mapping_from_groupid_to_set = self.mapping_from_groupid_to_set
    
        if not referenced_group.id in mapping_from_groupid_to_set:
            linked_set = self.load_from_group(referenced_group)
        else:
            linked_set = mapping_from_groupid_to_set[referenced_group.id]
            
        return linked_set
        
    def derefence(self, reference):
        return self.hdf5file[reference]
        
    def new_version(self, master_group):
        index = len(master_group)
        name = format(index + 1,"010d")
        return master_group.create_group(name)
        
    def info_group(self, ensure = True):
        return self.named_group(self.INFO_GROUP_NAME, ensure)
        
    def data_group(self, ensure = True):
        return self.named_group(self.DATA_GROUP_NAME, ensure)
        
    def named_group(self, name, ensure = True):
        if self.hdf5file.mode == 'r' or not ensure:
            if not name in self.hdf5file:
                return None
            else:
                return self.hdf5file[name]
        else:
            return self.hdf5file.require_group(name)
        
    def close(self):
        if not self.hdf5file is None:
            self.hdf5file.flush()
            self.hdf5file.close()
            self.hdf5file = None


class HDF5UnicodeAttribute(HDF5UnitlessAttribute):
    
    def __init__(self, name, dataset):
        HDF5UnitlessAttribute.__init__(self, name, dataset)
        
    def get_values(self, indices):
        if indices is None:
            encoded = self.dataset[:]
        elif len(self.get_shape()) > 1:
            encoded = self.dataset[:][indices]
        else:
            encoded = self.dataset[:][indices]
        return numpy.char.decode(encoded, 'UTF-32BE')

    def set_values(self, indices, values):
        self.dataset[indices] = numpy.char.encode(values, 'UTF-32LE')
    
    def get_value(self, index):
        return self.dataset[index].decode('UTF-32BE')



