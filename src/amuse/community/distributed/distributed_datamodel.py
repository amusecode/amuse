from amuse.datamodel import Particles, Particle


class Resources(Particles):
    
    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]
        index = self.get_all_indices_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Resource(keys, self, index, self._get_version())
    
    def __iter__(self):
        keys =  self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        version = self._get_version()

        for i in range(len(keys)):
            yield Resource(keys[i], self,  indices[i], version)
    
    add_resource = Particles.add_particle
    add_resources = Particles.add_particles
    remove_resource = Particles.remove_particle
    remove_resources = Particles.remove_particles


class Resource(Particle):
    
    def __init__(self, key = None, particles_set = None, set_index = -1, set_version = -1, **keyword_arguments):
        if particles_set is None:
            if key == None:
                particles_set = Resources(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = Resources(1, keys = [key])

        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        object.__setattr__(self, "_set_index", set_index)
        object.__setattr__(self, "_set_version", set_version)

        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)




class Pilots(Particles):
    
    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]
        index = self.get_all_indices_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Pilot(keys, self, index, self._get_version())
    
    def __iter__(self):
        keys =  self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        version = self._get_version()

        for i in range(len(keys)):
            yield Pilot(keys[i], self,  indices[i], version)
    
    add_pilot = Particles.add_particle
    add_pilots = Particles.add_particles
    remove_pilot = Particles.remove_particle
    remove_pilots = Particles.remove_particles


class Pilot(Particle):
    
    def __init__(self, key = None, particles_set = None, set_index = -1, set_version = -1, **keyword_arguments):
        if particles_set is None:
            if key == None:
                particles_set = Pilots(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = Pilots(1, keys = [key])

        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        object.__setattr__(self, "_set_index", set_index)
        object.__setattr__(self, "_set_version", set_version)

        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
            
class ScriptJobs(Particles):
    
    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]
        index = self.get_all_indices_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return ScriptJob(keys, self, index, self._get_version())
    
    def __iter__(self):
        keys =  self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        version = self._get_version()

        for i in range(len(keys)):
            yield ScriptJob(keys[i], self,  indices[i], version)
    
    submit_script_job = Particles.add_particle
    submit_script_jobs = Particles.add_particles
    cancel_script_job = Particles.remove_particle
    cancel_script_jobs = Particles.remove_particles
    remove_script_job = Particles.remove_particle
    remove_script_jobs = Particles.remove_particles

class ScriptJob(Particle):
    
    def __init__(self, key = None, particles_set = None, set_index = -1, set_version = -1, **keyword_arguments):
        if particles_set is None:
            if key == None:
                particles_set = ScriptJobs(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = ScriptJobs(1, keys = [key])

        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        object.__setattr__(self, "_set_index", set_index)
        object.__setattr__(self, "_set_version", set_version)

        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
class FunctionJobs(Particles):
    
    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]
        index = self.get_all_indices_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return FunctionJob(keys, self, index, self._get_version())
    
    def __iter__(self):
        keys =  self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        version = self._get_version()

        for i in range(len(keys)):
            yield FunctionJob(keys[i], self,  indices[i], version)
    
    submit_function_job = Particles.add_particle
    submit_function_jobs = Particles.add_particles
    cancel_function_job = Particles.remove_particle
    cancel_function_jobs = Particles.remove_particles
    remove_function_job = Particles.remove_particle
    remove_function_jobs = Particles.remove_particles

class FunctionJob(Particle):
    
    def __init__(self, key = None, particles_set = None, set_index = -1, set_version = -1, **keyword_arguments):
        if particles_set is None:
            if key == None:
                particles_set = FunctionJobs(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = FunctionJobs(1, keys = [key])

        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        object.__setattr__(self, "_set_index", set_index)
        object.__setattr__(self, "_set_version", set_version)

        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)


