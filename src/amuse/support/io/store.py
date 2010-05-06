import h5py
import numpy
import pickle
import os.path

from amuse.support.data.core import Particles, AttributeStorage
from amuse.support.units import si, units, core
from amuse.support.io import base


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
        
    def _get_attributes(self):
        return self.attributesgroup.keys()
        
    def _get_values(self, particles, attributes):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute in attributes:
             values_vector = self.attributesgroup[attribute]
             selection = h5py.selections.PointSelection(values_vector.shape)
             selection.set(numpy.transpose([indices,]))
             selected_values = values_vector[selection]
             results.append(self.get_unit_of(attribute).new_quantity(selected_values))
        
        return results
        
    def _has_key(self, key):
        return key in self.mapping_from_particle_to_index
    
    def _get_keys(self):
        return self.particle_keys
        
    def _set_values(self, particles, attributes, quantities):
        indices = self.get_indices_of(particles)
        
        for attribute, quantity in zip(attributes, quantities):
            if attribute in self.attributesgroup:
                dataset = self.attributesgroup[attribute]
            else:
                dataset = self.attributesgroup.create_dataset(attribute, shape=len(self.particle_keys), dtype=quantity.number.dtype)
                dataset["unit"] =  quantity.unit.to_simple_form().reference_string()
            dataset[indices] = quantity.value_in(self.get_unit_of(attribute))
        
        
class StoreHDF(object):
    PARTICLES_GROUP_NAME = "particles"
    
    def __init__(self, filename, append_to_file=True):
        if not append_to_file and os.path.exists(filename):
            os.remove(filename)
        self.hdf5file = h5py.File(filename ,'a')
        
    def store(self, particles):
        particles_group = self.particles_group()
        index = len(particles_group)
        name = format(index + 1,"010d")
        
        group = particles_group.create_group(name)
        
        attributes_group = group.create_group("attributes")
        
        
        group.attrs["number_of_particles"] = len(particles)
        group.attrs["class_of_the_particles"] = pickle.dumps(particles._particles_factory())
        
        quantity = particles.get_timestamp()
        if not quantity is None:
            group.attrs["timestamp"] = quantity.value_in(quantity.unit)
            group.attrs["timestamp_unit"] = quantity.unit.reference_string()
        
        #times_group = group.create_group("times")
        #previous_model_times = None
        #all_stored_model_time_references = {}
        #for attribute_values in particles._get_attributes():
        #    index = 0
        #    if not attribute_values.model_times is previous_model_times:
        #        model_times = attribute_values.model_times
        #        dataset = times_group.create_dataset(str(index), data=model_times.number)
        #        dataset.attrs["units"] = str(model_times.unit.to_simple_form())
        #        all_stored_model_time_references[id(model_times)] = index
        #        index += 1
        keys = particles._get_keys()
        dataset = group.create_dataset("keys", data=keys)
        
        all_values = particles._get_values(None, particles._get_attributes())
    
        for attribute, quantity in zip(particles._get_attributes(), all_values):
            value = quantity.value_in(quantity.unit)
            dataset = attributes_group.create_dataset(attribute, data=value)
            dataset.attrs["units"] = quantity.unit.to_simple_form().reference_string()
            #if not attribute_values.model_times is None:
            #    dataset.attrs["model-times-ref"] = all_stored_model_time_references[id(attribute_values.model_times)]
    
        self.hdf5file.flush()
        
    def load(self):
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
            
        return all_particle_sets[-1]
        
        
    
    def particles_group(self):
        return self.hdf5file.require_group(self.PARTICLES_GROUP_NAME)
        
    def close(self):
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
    
