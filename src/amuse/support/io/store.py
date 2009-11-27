import h5py

class StoreHDF(object):
    
    def __init__(self, filename):
        self.h5file = h5py.File(filename)
        
    def store(self, particles):
        group = self.h5file.create_group("amuse")
        
        attributes_group = group.create_group("attributes")
        times_group = group.create_group("times")
        
        previous_model_times = None
        all_stored_model_time_references = {}
        for attribute_values in particles.attributelist.mapping_from_attribute_to_values_and_unit.values():
            index = 0
            if not attribute_values.model_times is previous_model_times:
                model_times = attribute_values.model_times
                dataset = attributes_group.create_dataset(str(index), data=model_times.number)
                dataset.attrs["units-si"] = str(model_times.unit.to_simple_form())
                all_stored_model_time_references[id(model_times)] = index
                
        for attribute, attribute_values in particles.attributelist.mapping_from_attribute_to_values_and_unit.iteritems():
            dataset = attributes_group.create_dataset(attribute, data=attribute_values.values)
            dataset.attrs["units-si"] = str(attribute_values.unit.to_simple_form())
            if not attribute_values.model_times is None:
                dataset.attrs["model-times-ref"] = all_stored_model_time_references[id(attribute_values.model_times)]



