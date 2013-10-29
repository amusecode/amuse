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




class Reservations(Particles):
    
    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]
        index = self.get_all_indices_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Reservation(keys, self, index, self._get_version())
    
    def __iter__(self):
        keys =  self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        version = self._get_version()

        for i in range(len(keys)):
            yield Reservation(keys[i], self,  indices[i], version)
    
    add_reservation = Particles.add_particle
    add_reservations = Particles.add_particles


class Reservation(Particle):
    
    def __init__(self, key = None, particles_set = None, set_index = -1, set_version = -1, **keyword_arguments):
        if particles_set is None:
            if key == None:
                particles_set = Reservations(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = Reservations(1, keys = [key])

        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        object.__setattr__(self, "_set_index", set_index)
        object.__setattr__(self, "_set_version", set_version)

        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)


