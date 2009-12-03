import h5py
import numpy
import pickle
import os.path

from amuse.support.data.core import AttributeList, Particles
from amuse.support.units import si, units

class StoreHDF(object):
    PARTICLES_GROUP_NAME = "particles"
    
    def __init__(self, filename):
        if os.path.exists(filename):
            self.hdf5file = h5py.File(filename ,'r+')
        else:
            self.hdf5file = h5py.File(filename ,'w')
        
    def store(self, particles):
        particles_group = self.particles_group()
        index = len(particles_group)
        name = format(index + 1,"010d")
        
        group = particles_group.create_group(name)
        
        attributes_group = group.create_group("attributes")
        times_group = group.create_group("times")
        
        
        group.attrs["number_of_particles"] = len(particles)
        group.attrs["class_of_the_particles"] = pickle.dumps(type(particles))
        
        previous_model_times = None
        all_stored_model_time_references = {}
        for attribute_values in particles.attributelist.mapping_from_attribute_to_values_and_unit.values():
            index = 0
            if not attribute_values.model_times is previous_model_times:
                model_times = attribute_values.model_times
                dataset = times_group.create_dataset(str(index), data=model_times.number)
                dataset.attrs["units"] = str(model_times.unit.to_simple_form())
                all_stored_model_time_references[id(model_times)] = index
                index += 1
                
        for attribute, attribute_values in particles.attributelist.mapping_from_attribute_to_values_and_unit.iteritems():
            dataset = attributes_group.create_dataset(attribute, data=attribute_values.values)
            dataset.attrs["units"] = str(attribute_values.unit.to_simple_form())
            if not attribute_values.model_times is None:
                dataset.attrs["model-times-ref"] = all_stored_model_time_references[id(attribute_values.model_times)]
    
    def load(self):
        particles_group = self.particles_group()
        number_of_saved_particle_sets = len(particles_group)
        all_particle_sets = [None] * number_of_saved_particle_sets
        for group_index in particles_group.keys():
            group = particles_group[group_index]
            attributes_group = group["attributes"]
            times_group = group["times"]
            number_of_particles = group.attrs["number_of_particles"]
            class_of_the_particles = pickle.loads(group.attrs["class_of_the_particles"])
            particles = class_of_the_particles(number_of_particles)
            all_units = []
            lists_of_values = []
            attributes = []
            for attribute in attributes_group.keys():
                dataset = attributes_group[attribute]
                array = numpy.ndarray(len(dataset), dtype = dataset.dtype)
                all_units.append(eval(dataset.attrs["units"], units.__dict__))
                dataset.read_direct(array)
                lists_of_values.append(array)
                attributes.append(attribute)
            
            model_time_lists = [None] * len(times_group)
            for time_index in times_group.keys():
                dataset = times_group[time_index ]
                unit = eval(dataset.attrs["units"], units.__dict__)
                array = numpy.ndarray(len(dataset), dtype = dataset.dtype)
                dataset.read_direct(array)
                
                model_time = unit.new_quantity(array)
                model_time_lists[int(time_index)] = model_time
            
            
            particles.attributelist = AttributeList(
                    particles.particle_keys,
                    attributes = attributes, 
                    lists_of_values = lists_of_values, 
                    units = all_units)
                    
            for attribute, attribute_values in particles.attributelist.mapping_from_attribute_to_values_and_unit.iteritems():
                dataset = attributes_group[attribute]
                if "model-times-ref" in dataset.attrs:
                    index = int(dataset.attrs["model-times-ref"])
                    attribute_values.model_times = model_time_lists[index]
                    
                    
            all_particle_sets[int(group_index) - 1] = particles\
            
        previous = None
        for x in all_particle_sets:
            x.previous = previous
            previous = x
            
        return all_particle_sets[-1]
        
        
    
    def particles_group(self):
        return self.hdf5file.require_group(self.PARTICLES_GROUP_NAME)
            

