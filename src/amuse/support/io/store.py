import h5py
import numpy
import pickle
import os.path

from amuse.support.data.core import Particles, AttributeStorage
from amuse.support.units import si, units, core
from amuse.support.io import base

from amuse.support import exceptions


class HDF5AttributeStorage(AttributeStorage):

    def __init__(self, keys, hdfgroup):
        self.hdfgroup = hdfgroup
        self.attributesgroup = self.hdfgroup["attributes"]
        self.number_of_particles = self.hdfgroup.attrs["number_of_particles"]
        self.particle_keys = keys
        self.mapping_from_particle_to_index = self.new_index()
    
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
    
    def get_unit_of(self, attribute):
        dataset = self.attributesgroup[attribute]
        return eval(dataset.attrs["units"], core.__dict__) 
        
    def get_attribute_names_defined_in_store(self):
        return self.attributesgroup.keys()
        
    def get_values_in_store(self, particles, attributes):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute in attributes:
            values_vector = self.attributesgroup[attribute]
            selection = h5py.selections.PointSelection(values_vector.shape)
            selection.set(numpy.transpose([indices,]))
            selected_values = values_vector[selection]
            results.append(self.get_unit_of(attribute).new_quantity(selected_values))
        
        return results
        
    def has_key_in_store(self, key):
        return key in self.mapping_from_particle_to_index
    
    def get_all_keys_in_store(self):
        return self.particle_keys
        
    def set_values_in_store(self, particles, attributes, quantities):
        indices = self.get_indices_of(particles)
        
        for attribute, quantity in zip(attributes, quantities):
            if attribute in self.attributesgroup:
                dataset = self.attributesgroup[attribute]
            else:
                dataset = self.attributesgroup.create_dataset(attribute, shape=len(self.particle_keys), dtype=quantity.number.dtype)
                dataset["unit"] =  quantity.unit.to_simple_form().reference_string()
            dataset[indices] = quantity.value_in(self.get_unit_of(attribute))


class HDF5GridAttributeStorage(AttributeStorage):

    def __init__(self, shape, hdfgroup):
        self.hdfgroup = hdfgroup
        self.attributesgroup = self.hdfgroup["attributes"]
        self.shape = shape
        
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
        
    def get_attribute_names_defined_in_store(self):
        return self.attributesgroup.keys()
        
    def get_values_in_store(self, indices, attributes):
            
        results = []
        for attribute in attributes:
            values_vector = self.attributesgroup[attribute]
            if indices is None:
                selected_values = numpy.ndarray(shape=self.shape, dtype = values_vector.dtype)
                values_vector.read_direct(selected_values)
            else:    
                #selection = h5py.selections.PointSelection(values_vector.shape)
                #selection.set(numpy.transpose([indices,]))
                selected_values = values_vector[indices]
            results.append(self.get_unit_of(attribute).new_quantity(selected_values))
        
        return results
        
    def has_key_in_store(self, key):
        return False
    
    def get_all_keys_in_store(self):
        return None
        
    def set_values_in_store(self, indices, attributes, quantities):
        
        for attribute, quantity in zip(attributes, quantities):
            if attribute in self.attributesgroup:
                dataset = self.attributesgroup[attribute]
            else:
                dataset = self.attributesgroup.create_dataset(attribute, shape=self.shape, dtype=quantity.number.dtype)
                dataset["unit"] =  quantity.unit.to_simple_form().reference_string()
            dataset[indices] = quantity.value_in(self.get_unit_of(attribute))
        
        
class StoreHDF(object):
    PARTICLES_GROUP_NAME = "particles"
    GRIDS_GROUP_NAME = "grids"
    
    def __init__(self, filename, append_to_file=True):
        if not append_to_file and os.path.exists(filename):
            os.remove(filename)
        self.hdf5file = h5py.File(filename,'a')
    
    
    def store(self, container):
        if hasattr(container, 'shape') and len(container.shape) > 1:
            self.store_grid(container)
        else:            
            self.store_particles(container)
    
    def new_group(self, master_group):
        index = len(master_group)
        name = format(index + 1,"010d")
        return master_group.create_group(name)
        
    def store_particles(self, particles):
        group = self.new_group(self.particles_group())
        
        group.attrs["number_of_particles"] = len(particles)
        group.attrs["class_of_the_particles"] = pickle.dumps(particles._factory_for_new_collection())
            
        keys = particles.get_all_keys_in_store()
        dataset = group.create_dataset("keys", data=keys)

        self.store_timestamp(particles, group)
        self.store_values(particles, group)
            
        self.hdf5file.flush()
    
    def store_grid(self, grid):
        group = self.new_group(self.grids_group())
        
        group.attrs["class_of_the_container"] = pickle.dumps(grid._factory_for_new_collection())
        group.create_dataset("shape", data=numpy.asarray(grid.shape))
        
        self.store_timestamp(grid, group)
        self.store_values(grid, group)
        
        self.hdf5file.flush()
        
    def store_values(self, container, group):
        attributes_group = group.create_group("attributes")
        
        all_values = container.get_values_in_store(None, container.get_attribute_names_defined_in_store())
        for attribute, quantity in zip(container.get_attribute_names_defined_in_store(), all_values):
            value = quantity.value_in(quantity.unit)
            dataset = attributes_group.create_dataset(attribute, data=value)
            dataset.attrs["units"] = quantity.unit.to_simple_form().reference_string()
            
    
    def store_timestamp(self, container, group):
        quantity = container.get_timestamp()
        if not quantity is None:
            group.attrs["timestamp"] = quantity.value_in(quantity.unit)
            group.attrs["timestamp_unit"] = quantity.unit.reference_string()
    
    
    def load(self):
        if not self.PARTICLES_GROUP_NAME in self.hdf5file:
            return self.load_grid()
        else:
            if len(self.particles_group()) > 0:
                return self.load_particles()
            else:
                return self.load_grid()
                
    def load_particles(self):
        particles_group = self.particles_group()
        number_of_saved_particle_sets = len(particles_group)
        all_particle_sets = [None] * number_of_saved_particle_sets
        for group_index in particles_group.keys():
            group = particles_group[group_index]
            class_of_the_particles = pickle.loads(group.attrs["class_of_the_particles"])
            dataset = group["keys"]
            keys = numpy.ndarray(len(dataset), dtype = dataset.dtype)
            dataset.read_direct(keys)
            
            if "timestamp" in group.attrs:
                unit = eval(group.attrs["timestamp_unit"], core.__dict__) 
                timestamp = unit.new_quantity(group.attrs["timestamp"])
            else:
                timestamp = None
            
            particles = class_of_the_particles()
            particles._private.attribute_storage = HDF5AttributeStorage(keys, group)
            particles._private.timestamp = timestamp
            
            all_particle_sets[int(group_index) - 1] = particles
            
        previous = None
        for x in all_particle_sets:
            x._private.previous = previous
            previous = x
            
        last = all_particle_sets[-1]
        copy_of_last = last.copy_to_memory()
        copy_of_last._private.previous = last
        return copy_of_last
        
        
    
    def load_grid(self):
        container_group = self.grids_group()
        number_of_saved_particle_containers= len(container_group)
        all_containers = [None] * number_of_saved_particle_containers
        for group_index in container_group.keys():
            group = container_group[group_index]
            class_of_the_container = pickle.loads(group.attrs["class_of_the_container"])
            shape = tuple(group["shape"])
            
            if "timestamp" in group.attrs:
                unit = eval(group.attrs["timestamp_unit"], core.__dict__) 
                timestamp = unit.new_quantity(group.attrs["timestamp"])
            else:
                timestamp = None
            
            container = class_of_the_container()
            container._private.attribute_storage = HDF5GridAttributeStorage(shape, group)
            container._private.timestamp = timestamp
            
            all_containers[int(group_index) - 1] = container
            
        previous = None
        for x in all_containers:
            x._private.previous = previous
            previous = x
            
        last = all_containers[-1]
        copy_of_last = last.copy_to_memory()
        copy_of_last._private.previous = last
        return copy_of_last
        
    def particles_group(self):
        return self.hdf5file.require_group(self.PARTICLES_GROUP_NAME)
        
    def grids_group(self):
        return self.hdf5file.require_group(self.GRIDS_GROUP_NAME)
        
    def close(self):
        self.hdf5file.flush()
        self.hdf5file.close()

class HDF5FileFormatProcessor(base.FileFormatProcessor):
    """
    Process an HDF5 file
    """
    
    provided_formats = ['hdf5', 'amuse']
    
    def __init__(self, filename = None, set = None, format = None, append_to_file=True):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        self.append_to_file = append_to_file
    
    def load(self):
        processor = StoreHDF(self.filename, self.append_to_file)
        return processor.load()
        
    def store(self):
        processor = StoreHDF(self.filename, self.append_to_file)
        return processor.store(self.set)
    
    @base.format_option
    def append_to_file(self):
        "By default new data is appended to HDF5 files. Set this to False to overwrite existing files."
        return True
    
