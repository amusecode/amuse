"""
"""

from amuse.support.data import values
from amuse.support.units import si
from amuse.support.units import units

import numpy
import random

class BasicUniqueKeyGenerator(object):
    
    def __init__(self):
        self.lowest_unique_key = 1
    
    def next(self):
        new_key = self.lowest_unique_key
        self.lowest_unique_key += 1
        return new_key
        
    def next_set_of_keys(self, length):
        if length == 0:
            return  []
            
        from_key = self.lowest_unique_key
        to_key = from_key + length;
        self.lowest_unique_key += length
        return numpy.arange(from_key, to_key)
        

class RandomNumberUniqueKeyGenerator(object):
    DEFAULT_NUMBER_OF_BITS = 64
    
    def __init__(self, number_of_bits = DEFAULT_NUMBER_OF_BITS):
        if number_of_bits > 64:
            raise Exception("number of bits is larger than 64, this is currently unsupported!")
        self.number_of_bits = number_of_bits
        
    def next(self):
        return random.getrandbits(self.number_of_bits)
        
    def next_set_of_keys(self, length):
        if length == 0:
            return  []
        return numpy.array([random.getrandbits(64) for i in range(length)], dtype='uint64')
        
UniqueKeyGenerator = RandomNumberUniqueKeyGenerator()

class AttributeValues(object):
    __slots__ = ["attribute", "values", "unit", "model_times"]
    
    def __init__(self, attribute, unit, values = None,  model_times = None, length = None):
        self.attribute = attribute
        self.unit = unit
        self.model_times = model_times
        if values is None:
            self.values = numpy.zeros(length)
        else:
            self.values = values
        
    def copy(self):
        return AttributeValues(self.attribute, self.unit, self.values.copy(), self.model_times)
        
class InMemoryAttributeStorage(object):
    
    def __init__(self):
        self.model_times = None
        self.mapping_from_attribute_to_values_and_unit = {}
        self.mapping_from_particle_to_index = {}
        self.particle_keys = []

    def _set_particles(self, keys, attributes = [], values = []):
        if len(values) != len(attributes):
            raise Exception(
                "you need to provide the same number of value list as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(values)
                )
            )
        if len(values) > 0 and len(keys) != len(values[0]):
            raise Exception(
                "you need to provide the same number of values as particles, found {0} values and {1} particles".format(
                    len(values[0]), len(keys)
                )
            )
        
        if len(self.particle_keys) > 0:
            for attribute, values_to_set in zip(attributes, values):
                if attribute in self.mapping_from_attribute_to_values_and_unit:
                    attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
                else:
                    attribute_values = AttributeValues(
                        attribute,
                        values_to_set.unit,
                        length = len(self.particle_keys)
                    )
                
                    self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                values_in_the_right_units = values_to_set.value_in(attribute_values.unit)
                attribute_values.values = numpy.concatenate((attribute_values.values, values_in_the_right_units))
            
            old_length = len(self.particle_keys)
            zeros_for_concatenation = numpy.zeros(len(keys))
            for attribute_values in self.mapping_from_attribute_to_values_and_unit.values():
                if len(attribute_values.values) == old_length:
                    attribute_values.values = numpy.concatenate((attribute_values.values, zeros_for_concatenation))
            
                    
            index = len(self.particle_keys)

            self.particle_keys = numpy.concatenate((self.particle_keys,  numpy.array(list(keys))))
  
            for particle_key in keys:
                self.mapping_from_particle_to_index[particle_key] = index
                index += 1
        else:
            self.mapping_from_attribute_to_values_and_unit = {}
            
            for attribute, quantity in zip(attributes, values):
                self.mapping_from_attribute_to_values_and_unit[attribute] = AttributeValues(
                    attribute,
                    quantity.unit,
                    quantity.number,
                    None
                )
             
            self.particle_keys = numpy.array(keys)
            
            self.reindex()
    
    def _get_values(self, particles, attributes):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute in attributes:
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             if indices is None:
                 selected_values = attribute_values.values
             else:
                 selected_values = attribute_values.values.take(indices)
             results.append(attribute_values.unit.new_quantity(selected_values))
        
        return results
        
    def _set_values(self, particles, attributes, list_of_values_to_set, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, len(indices))
        
        previous_model_times = None
        if list_of_values_to_set is None:
            for attribute in attributes:
                if attribute in self.mapping_from_attribute_to_values_and_unit:
                    attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
                else:
                    raise Exception("unknown attribute '{0}'".format(attribute))
                     
                selected_values = numpy.zeros(len(indices))
                
                attribute_values.values.put(indices, selected_values)
            return
            
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                   attribute,
                   values_to_set.unit,
                   length = len(self.particle_keys)
                )
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                 
            selected_values = values_to_set.value_in(attribute_values.unit)
            
            attribute_values.values.put(indices, selected_values)
            if not model_times is None:
                if not previous_model_times is attribute_values.model_times:
                    attribute_values.model_times.put(indices, model_times)
                    previous_model_times = attribute_values.model_times
            
    
    
    def _get_attributes(self):
        return sorted(self.mapping_from_attribute_to_values_and_unit.keys())
    
    
    def _has_key(self, key):
        return key in self.mapping_from_particle_to_index
        
    def _get_keys(self):
        return self.particle_keys
        
    def __len__(self):
        return len(self.particle_keys)
        
    def copy(self):
        copy = InMemoryAttributeStorage()
        copy.mapping_from_particle_to_index = self.mapping_from_particle_to_index.copy()
        copy.particle_keys = self.particle_keys.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_values_and_unit.iteritems():
            copy.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle_key, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            raise AttributeError("particle does not have a "+attribute)
        
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        index = self.mapping_from_particle_to_index[particle_key]
        
        return attribute_values.unit.new_quantity(attribute_values.values[index])
        
    def set_value_of(self, particle_key, attribute, quantity):
        index = self.mapping_from_particle_to_index[particle_key]
            
        if index is None:
            raise Exception("particle with key '{0}' is not in this set".format(particle_key))
            
        attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
             
        value_to_set = quantity.value_in(attribute_values.unit)
        attribute_values.values[index] = value_to_set
        
            
    def iter_values_of_particle(self, particle_key):
        index = self.mapping_from_particle_to_index[particle_key]
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            yield attribute, (attribute_values.values[index] | attribute_values.unit)
    
    
            
    def iter_values_of(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            return
            
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        values = attribute_values.values
        unit = attribute_values.unit
        particles = self.particle_keys
        
        for index in range(len(self.particle_keys)):
            yield particles[i], (values[i] | unit)
            
   
    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
            
        mapping_from_particle_to_index = self.mapping_from_particle_to_index 
        result = numpy.zeros(len(particles),dtype='int32')
        #result = [mapping_from_particle_to_index[id(particle)] for particle in particles]
        
        index = 0
        for index, particle_key in enumerate(particles):
            result[index] = mapping_from_particle_to_index[particle_key]
            index += 1
        return result
    
    def get_values_of_particles_in_units(self, particles, attributes, target_units):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute, target_unit in zip(attributes, target_units):
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_unit_in_target_unit = attribute_values.unit.value_in(target_unit )
             if indices is None:
                 selected_values = attribute_values.values
             else:
                 selected_values = attribute_values.values.take(indices)
             if value_of_unit_in_target_unit != 1.0:
                selected_values *= value_of_unit_in_target_unit
             results.append(selected_values)
        
        return results
        
    
            
    def _convert_model_times(self, model_times, length):
        if not model_times is None and isinstance(model_times, values.ScalarQuantity):
            return model_times.unit.new_quantity(numpy.linspace(model_times.number, model_times.number, length) )
        else:
            return model_times
    
    def set_values_of_particles_in_units(self, particles, attributes, list_of_values_to_set, source_units, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, len(indices))
        
        previous_model_times = None
        results = []
        for attribute, values_to_set, source_unit in zip(attributes, list_of_values_to_set, source_units):
            selected_values = values_to_set
            
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                   attribute,
                   source_unit,
                   length = len(self.particle_keys)
                )
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                 
            value_of_source_unit_in_list_unit = source_unit.value_in(attribute_values.unit)
            if value_of_source_unit_in_list_unit != 1.0:
                selected_values *= value_of_source_unit_in_list_unit 
             
            attribute_values.values.put(indices, selected_values)
            if not model_times is None:
                if not previous_model_times is attribute_values.model_times:
                    attribute_values.model_times.put(indices, model_times)
                    previous_model_times = attribute_values.model_times
            
        return results
        
             
    def merge_into(self, others):
        source_attributes = []
        source_units = []
        source_valeus = []
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            source_attributes.append(attribute_values.attribute)
            source_values.append(attribute_values.values)
            source_units.append(attribute_values.unit)
            
                
        other.set_values_of_particles_in_units(self.particle_keys, source_attributes, source_values, source_units)
        
    def remove_particles(self, particles):
        indices = self.get_indices_of(particles)
        
        mapping_from_attribute_to_values_and_unit = self.mapping_from_attribute_to_values_and_unit.copy()
        for attribute, attribute_values in mapping_from_attribute_to_values_and_unit.iteritems():
            attribute_values.values = numpy.delete(attribute_values.values,indices)
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
        
    def _remove_particles(self, keys):
        indices = self.get_indices_of(keys)
        
        for attribute, attribute_values in self.mapping_from_attribute_to_values_and_unit.iteritems():
            attribute_values.values = numpy.delete(attribute_values.values,indices)
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
        
    def reindex(self):
        d = {}
        index = 0
        for particle_key in self.particle_keys:
            d[particle_key] = index
            index += 1
          
        self.mapping_from_particle_to_index = d
        
    def get_attribute_as_vector(self, particle_key, names_of_the_scalar_attributes):
        index = self.mapping_from_particle_to_index[particle_key]
        vector = []
        for i, attribute in enumerate(names_of_the_scalar_attributes):
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            vector.append(attribute_values.values[index])
            unit_of_the_vector = attribute_values.unit
        return vector | unit_of_the_vector
                
    def set_attribute_as_vector(self, particle_key, names_of_the_scalar_attributes, quantity):
        index = self.mapping_from_particle_to_index[particle_key]
        vector = []
        for i, attribute  in enumerate(names_of_the_scalar_attributes):
            attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
        
            value_to_set = quantity[i].value_in(attribute_values.unit)
            attribute_values.values[index] = value_to_set
           
            
    
    def get_attributevalues_for(self, attribute, unit):
        if attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        else:
            attribute_values = AttributeValues(
                attribute,
                unit,
                length = len(self.particle_keys)
            )
            self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
        return attribute_values


    def attributes(self):
        return set(self.mapping_from_attribute_to_values_and_unit.keys())
    
    def __str__(self):
        attributes = sorted(self.attributes())
        
        columns = map(lambda x : [x], attributes)
        columns.insert(0,['id'])
        
        for index, attribute in enumerate(attributes):
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            column = columns[index + 1]
            column.append(str(attribute_values.unit))
            column.append('========')
            if len(attribute_values.values) > 40:
                values_to_show = list(attribute_values.values[:20])
                values_to_show.append(attribute_values.values[-20:])
            else:
                values_to_show = attribute_values.values
            
            for value in values_to_show:
                column.append(str(value))
            column.append('========')
            
        column = columns[0]
        column.append("-")
        column.append('========')
        if len(self.particle_keys) > 40:
            values_to_show = list(self.particle_keys[:20])
            values_to_show.append(self.particle_keys[-20:])
        else:
            values_to_show = self.particle_keys
                    
        for value in values_to_show:
            column.append(str(value))
            
        column.append('========')
            
        rows = []
        for i in range(len(columns[0])):
            row = [x[i] for x in columns]
            rows.append(row)
            
        lines = map(lambda  x : '\t'.join(x), rows)
        return '\n'.join(lines)        
        
        
        
    def get_attributes_as_vector_quantities(self, attributes):
        for attribute in attributes:
            if not attribute in self.mapping_from_attribute_to_values_and_unit:
                raise AttributeError("unknown attribute '{0}'".format(attribute))
        
        unit_of_the_values = None
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            if unit_of_the_values is None:
                unit_of_the_values = attribute_values.unit
                results.append(attribute_values.values)
            else:
                conversion_factor = attribute_values.unit.value_in(unit_of_the_values)
                if conversion_factor != 1.0:
                    results.append(attribute_values.values * conversion_factor)
                else:                    
                    results.append(attribute_values.values)
        
        results = numpy.dstack(results)[0]
        return unit_of_the_values.new_quantity(results)
        
    def set_model_time(self, value): 
        model_times = self._convert_model_times(value, len(self.particle_keys))
        for attribute_values in self.mapping_from_attribute_to_values_and_unit.values():
            attribute_values.model_times = model_times


class AbstractParticleSet(object):
    
    class PrivateProperties(object):
        """
        Defined for superclasses to store private properties.
        As __setattr__ is defined on subclasses these
        do not need to use object.__setattr__ syntax.
        """
        pass
    
    class VectorAttribute:
        
        def  __init__(self, attribute_names):
            self.attribute_names = attribute_names
        
        def _get_values(self, instance):
            return instance._private.attribute_storage.get_attributes_as_vector_quantities(self.attribute_names)
                
        def _set_values(self, instance, value):
            vectors = value.number
            split = numpy.hsplit(vectors,len(self.attribute_names))
            list_of_values = []
            for i in range(len(self.attribute_names)):
                values = value.unit.new_quantity(split[i].reshape(len(split[i])))
                list_of_values.append(values)
                
            instance._set_values(instance._get_keys(), self.attribute_names, list_of_values)
    
    
        
    def __init__(self):
        object.__setattr__(self, "_vector_attributes", {})
        object.__setattr__(self, "_private", self.PrivateProperties())
    
    def __getattr__(self, name_of_the_attribute):
        if name_of_the_attribute in self._vector_attributes:
            return self._vector_attributes[name_of_the_attribute]._get_values(self)
        else:
            return self._get_values(self._get_keys(), [name_of_the_attribute])[0]
    
    def __setattr__(self, name_of_the_attribute, value):
        if name_of_the_attribute in self._vector_attributes:
            self._vector_attributes[name_of_the_attribute]._set_values(self, value)
        else:
            self._set_values(self._get_keys(), [name_of_the_attribute], [value])
    
    #
    # Particle storage interface
    #
    
    def _get_values(self, keys, attribues):
        pass
    
    def _set_values(self, keys, attributes, values):
        pass
        
    def _set_particles(self, keys, attributes, values):
        pass
        
    def _remove_particles(self, keys):
        pass
    
    def _get_keys(self):
        return []
        
    def _has_key(self):
        return False
        
    def add_vector_attribute(self, name_of_the_attribute, name_of_the_components):
        self._vector_attributes[name_of_the_attribute] = self.VectorAttribute(name_of_the_components)
        
        
    #
    # public API
    #
    def __iter__(self):
        index = 0
        for key in self._get_keys():
            p = Particle(key, self)
            yield p
            index += 1

    def __getitem__(self, index):
        return Particle(self._get_keys()[index], self)
        
    def __len__(self):
        return len(self._get_keys())
        
     
    def copy(self):
         attributes = self._get_attributes()
         keys = self._get_keys()
         values = self._get_values(None, attributes)
         result = Particles()
         result._set_particles(keys, attributes, values)
         return result
         
    
    def copy_values_of_attribute_to(self, attribute_name, particles):
        channel = self.new_channel_to(particles)
        channel.copy_attributes([attribute_name])
    
    def new_channel_to(self, other):
        return ParticleInformationChannel(self, other)
        
    def add_particles(self, particles):
        attributes = self._get_attributes()
        keys = particles._get_keys()
        values = particles._get_values(None, attributes)
        self._set_particles(keys, attributes, values)
    
    def add_particle(self, particle):
        self.add_particles(particle.as_set())
        
    def synchronize_to(self, other_particles):
        other_keys = set(other_particles._get_keys())
        my_keys = set(self._get_keys())
        added_keys = my_keys - other_keys
        removed_keys = other_keys - my_keys
        
        added_keys = list(added_keys)
        if added_keys:
            attributes = self._get_attributes()
            values = self._get_values(added_keys, attributes)
            other_particles._set_particles(added_keys, attributes, values)
        
        removed_keys = list(removed_keys)
        if removed_keys:
            other_particles._remove_particles(removed_keys)
        
    def copy_values_of_state_attributes_to(self, particles):
        channel = self.new_channel_to(particles)
        channel.copy_attributes(self._private.attribute_storage._state_attributes())   

    
class Particles(AbstractParticleSet):
    
    """A set of particle objects"""
    def __init__(self, size = 0):
        AbstractParticleSet.__init__(self)
        particle_keys = UniqueKeyGenerator.next_set_of_keys(size)
        self._private.attribute_storage = InMemoryAttributeStorage()
        self._set_particles(particle_keys)
        self._private.previous = None
        
    def savepoint(self):
        instance = type(self)()
        instance._private.attribute_storage = self._private.attribute_storage.copy()
        instance._private.previous = self._private.previous
        self._private._previous = instance
        
    def iter_history(self):
        current = self
        while not current is None:
            yield current
            current = current._private.previous
    
    @property
    def history(self):
        return reversed(list(self.iter_history()))
        
    def get_timeline_of_attribute(self, particle_key, attribute):
        timeline = []
        for x in self.history:
            timeline.append((None, x._private.attribute_storage.get_value_of(particle_key, attribute)))
        return timeline
                    
    def copy(self):
         attributes = self._private.attribute_storage._state_attributes()
         keys = self._get_keys()
         values = self._get_values(None, attributes)
         result = Particles()
         result._set_particles(keys, attributes, values)
         return result
        
    def _set_particles(self, keys, attributes = [], values = []):
        self._private.attribute_storage._set_particles(keys, attributes, values)
        
    def _remove_particles(self, keys):
        self._private.attribute_storage._remove_particles(keys)
    
    def _get_values(self, keys, attributes):
        return self._private.attribute_storage._get_values(keys, attributes)
        
    def _set_values(self, keys, attributes, values):
        self._private.attribute_storage._set_values(keys, attributes, values)
    
    def _get_attributes(self):
        return self._private.attribute_storage._get_attributes()
        
    def _get_keys(self):
        return self._private.attribute_storage._get_keys()
        
    def _has_key(self, key):
        return self._private.attribute_storage._has_key(key)
        
              
              

class ParticlesSubset(AbstractParticleSet):
    """A sub set of particle objects"""
    def __init__(self, particles, keys):
        AbstractParticleSet.__init__(self)
        
        self._private.particles = particles
        self._private.keys = numpy.array(keys)
        self._private.set_of_keys = set(keys)
              
        
    def _set_particles(self, keys, attributes = [], values = []):
        self._private.keys = numpy.concatenate((self.keys,  numpy.array(keys)))
        self._private.set_of_keys += set(keys)
        self._private.particles._set_particles(keys, attributes, values)
        
    def _remove_particles(self, keys):
        raise Exception("Cannot remove particles from a subset")
    
    def _get_values(self, keys, attributes):
        return self._private.particles._get_values(keys, attributes)
        
    def _set_values(self, keys, attributes, values):
        self._private.particles._set_values(keys, attributes, values)
    
    def _get_attributes(self):
        return self._private.particles._get_attributes()
        
    def _get_keys(self):
        return self._private.keys
        
    def _has_key(self, key):
        return key in self._private.set_of_keys
        
        


class ParticleInformationChannel(object):
    
    def __init__(self, from_particles, to_particles):
        self.from_particles = from_particles
        self.to_particles = to_particles
        self._reindex()
        
    def _reindex(self):
        self.keys = self.intersecting_keys()
        #speed-up:
        #self.from_indices = self.from_particles._get_indices_of(self.keys)
        #self.to_indices = self.to_particles._get_indices_of(self.keys)
    
    def intersecting_keys(self):
        from_keys = self.from_particles._get_keys()
        return filter(lambda x : self.to_particles._has_key(x), from_keys)
        
    def copy_attributes(self, attributes):
        data = self.from_particles._get_values(self.keys, attributes)
        self.to_particles._set_values(self.keys, attributes, data)
    
class Stars(Particles):

    def __init__(self, size = 0):
        Particles.__init__(self, size)
        self.add_vector_attribute("position", ["x","y","z"])
        self.add_vector_attribute("velocity", ["vx","vy","vz"])
    
    
    def center_of_mass(self):
        masses = self.mass
        x_values = self.x
        y_values = self.y
        z_values = self.z
        
        total_mass = masses.sum()
        massx = (masses * x_values).sum()
        massy = (masses * y_values).sum()
        massz = (masses * z_values).sum()

        return values.VectorQuantity.new_from_scalar_quantities(
            massx/total_mass,
            massy/total_mass,
            massz/total_mass
        )

    
    
    def center_of_mass_velocity(self):
        masses = self.mass
        x_values = self.vx
        y_values = self.vy
        z_values = self.vz
        
        total_mass = masses.sum()
        massx = (masses * x_values).sum()
        massy = (masses * y_values).sum()
        massz = (masses * z_values).sum()

        return values.VectorQuantity.new_from_scalar_quantities(
            massx/total_mass,
            massy/total_mass,
            massz/total_mass
        )
        
    def kinetic_energy(self):
        mass = self.mass
        vx = self.vx
        vy = self.vy
        vz = self.vz
        v_squared = (vx * vx) + (vy * vy) + (vz * vz)
        m_v_squared = mass * v_squared
        return 0.5 * m_v_squared.sum()
        
    
    def potential_energy(self, smoothing_length_squared = 0 | units.m * units.m):
        mass = self.mass
        x_vector = self.x
        y_vector = self.y
        z_vector = self.z
        
        sum_of_energies = 0 | units.J
        
        for i in range(len(self)):
           x = x_vector[i]
           y = y_vector[i]
           z = z_vector[i]
           dx = x - x_vector[i+1:]
           dy = y - y_vector[i+1:]
           dz = z - z_vector[i+1:]
           dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
           dr = (dr_squared+smoothing_length_squared).sqrt()
           m_m = mass[i] * mass[i+1:]
           
           energy_of_this_particle = (m_m / dr).sum()
           sum_of_energies +=  -1 * units.G * energy_of_this_particle
            
        return sum_of_energies
        
        
        
        
        
            
        
            
class VectorAttribute(object):
    def __init__(self, names_of_the_scalar_components):
        self.names_of_the_scalar_components = names_of_the_scalar_components
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            return instance.particles_set._private.attribute_storage.get_attribute_as_vector(instance.key, self.names_of_the_scalar_components)
    
    
    def __set__(self, instance, quantity):
        instance.particles_set._private.attribute_storage.set_attribute_as_vector(instance.key, self.names_of_the_scalar_components, quantity)
            
class Particle(object):
    """A physical object or a physical region simulated as a 
    physical object (cload particle).
    
    All attributes defined on a particle are specific for 
    that particle (for example mass or position). A particle contains 
    a set of attributes, some attributes are *generic* and applicaple
    for multiple modules. Other attributes are *specific* and are 
    only applicable for a single module.
    """
    __slots__ = ["key", "particles_set"]
    
    
    position = VectorAttribute(("x","y","z"))
    velocity = VectorAttribute(("vx","vy","vz"))
    
    def __init__(self, key = None, particles_set = None, **keyword_arguments):
        if particles_set is None:
            particles_set = Particles(1)
            key = particles_set._get_keys()[0]
        
        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if hasattr(type(self), name_of_the_attribute):
            getattr(type(self), name_of_the_attribute).__set__(self, new_value_for_the_attribute)
            return
            
        if isinstance(new_value_for_the_attribute, values.Quantity):
            self.particles_set._private.attribute_storage.set_value_of(self.key, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise Exception("attribute "+name_of_the_attribute+" does not have a valid value, values must have a unit")
    
    def __getattr__(self, name_of_the_attribute):
         return self.particles_set._private.attribute_storage._get_values([self.key], [name_of_the_attribute])[0][0]
         
                
    def __str__(self):
        output = 'Particle '
        output += str(self.key)
        output += ''
        output += '\n'
        for name, value in self.particles_set._private.attribute_storage.iter_values_of_particle(self.key):
            output += name
            output += ': {'
            output += str(value)
            output += '}, '
            output += '\n'
        return output
        
    def set_default(self, attribute, quantity):
        if not attribute in self.set._private.attribute_storage.attributes():
            self.particles_set._private.attribute_storage.set_value_of(self, attribute, quantity)
            
    def get_timeline_of_attribute(self, attribute):
        return self.particles_set.get_timeline_of_attribute(self.key, attribute)
        
    def as_set(self):
        return ParticlesSubset(self.particles_set, [self.key])
            
