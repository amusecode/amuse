from amuse.support.data import values
from amuse.support.units import si
from amuse.support.units import units

import numpy
import random
import inspect

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


class AttributeStorage(object):
    
    def _set_particles(self, keys, attributes = [], values = []):
        pass
        
    def _remove_particles(self, keys):
        pass
        
    def _get_values(self, particles, attributes):
        pass
        
    def _set_values(self, particles, attributes, list_of_values_to_set):
        pass
        
    def _set_particles(self, keys, attributes = [], values = []):
        pass
        
    def _get_attributes(self):
        pass
    
    def _has_key(self, key):
        return False
        
    def _get_keys(self):
        return []
        
    def __len__(self):
        return 0
        
class InMemoryAttributeStorage(AttributeStorage):
    
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

    def attributes(self):
        return set(self.mapping_from_attribute_to_values_and_unit.keys())
    
    def _state_attributes(self):
        return self.attributes()
        
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
        
        
        
        
    def set_model_time(self, value): 
        model_times = self._convert_model_times(value, len(self.particle_keys))
        for attribute_values in self.mapping_from_attribute_to_values_and_unit.values():
            attribute_values.model_times = model_times




class AbstractParticleSet(object):
    """
    Abstract superclass of all sets of particles. 
    This class defines common code for all particle sets.
    
    Particle sets define dynamic attributes. Attributes
    can be set and retrieved on the particles using common python
    syntax. These attributes can only have values with units.
    
    >>> particles = Particles(2)
    >>> particles.mass = [10.0, 20.0] | units.kg
    >>> particles.mass
    quantity<[10.0, 20.0] kg>
    
    >>> particles.mass = 1.0 | units.kg
    >>> particles.mass
    quantity<[1.0, 1.0] kg>
    
    Particle sets can be iterated over. 
    
    >>> particles = Particles(2)
    >>> particles.mass = [10.0, 20.0] | units.kg
    >>> for particle in particles:
    ...     print particle.mass
    ...
    10.0 kg
    20.0 kg
    
    
    Particle sets can be indexed.
    
    >>> particles = Particles(3)
    >>> particles.x = [10.0, 20.0, 30.0] | units.m
    >>> particles[1].x
    quantity<20.0 m>
    
    
    Particle sets can be copied.
    
    >>> particles = Particles(3)
    >>> particles.x = [10.0, 20.0, 30.0] | units.m
    >>> copy = particles.copy()
    >>> particles.x = 2.0 | units.m
    >>> particles.x
    quantity<[2.0, 2.0, 2.0] m>
    >>> copy.x
    quantity<[10.0, 20.0, 30.0] m>
    
    Particle sets can have instance based or global vector attributes.
    A particle set stores a list of scalar values for each attribute.
    Some attributes are more naturally accessed as lists 
    of vector values. Once defined, a particle set can
    convert the scalar values of 2 or more attributes into one
    vector attribute.
    
    >>> particles = Particles(2)
    >>> particles.x = [1.0 , 2.0] | units.m
    >>> particles.y = [3.0 , 4.0] | units.m
    >>> particles.z = [5.0 , 6.0] | units.m
    >>> particles.add_vector_attribute("p", ["x","y","z"])
    >>> particles.p
    quantity<[[ 1.  3.  5.], [ 2.  4.  6.]] m>
    >>> particles.p[0]
    quantity<[1.0, 3.0, 5.0] m>
    >>> particles.position # "position" is a global vector attribute, coupled to x,y,z
    quantity<[[ 1.  3.  5.], [ 2.  4.  6.]] m>
    
    


    """
    GLOBAL_VECTOR_ATTRIBUTES = {}
    GLOBAL_CALCULATED_ATTRIBUTES = {}
    
    class PrivateProperties(object):
        """
        Defined for superclasses to store private properties.
        A particle set has ```__setattr__``` defined.
        The ```__setattr__``` function will set all attributes
        of the particles in the set to the specified value(s).
        To be able to define attributes on the set itself we
        use an instance of this class, attributes can be defined as::
        
            self._private.new_attribute = 'new value'
        
        Subclass implementers do not need to
        use the ```object.__setattr__``` syntax.
        
        For documentation about the __setattr__ call please
        see the ```python data model``` documentation on the python
        website.
        """
        pass
    
    class VectorAttribute:
        
        def  __init__(self, attribute_names):
            self.attribute_names = attribute_names
        
        def _get_values(self, instance):
            values = instance._get_values(instance._get_keys(), self.attribute_names)
              
            unit_of_the_values = None
            results = []
            for quantity in values:
                if unit_of_the_values is None:
                    unit_of_the_values = quantity.unit
                results.append(quantity.value_in(unit_of_the_values))
                
            results = numpy.dstack(results)[0]
            return unit_of_the_values.new_quantity(results)

        def _set_values(self, instance, value):
            vectors = value.number
            split = numpy.hsplit(vectors,len(self.attribute_names))
            list_of_values = []
            for i in range(len(self.attribute_names)):
                values = value.unit.new_quantity(split[i].reshape(len(split[i])))
                list_of_values.append(values)
                
            instance._set_values(instance._get_keys(), self.attribute_names, list_of_values)
        
        def get_value_for_particle(self, instance,  key):
            values = instance._get_values([key], self.attribute_names)
              
            unit_of_the_values = None
            results = []
            for quantity in values:
                if unit_of_the_values is None:
                    unit_of_the_values = quantity.unit
                results.append(quantity.value_in(unit_of_the_values))
                
            results = numpy.dstack(results)[0]
            return unit_of_the_values.new_quantity(results[0])
    
        def set_value_for_particle(self, instance, key, vector):
            list_of_values = []
            for quantity in vector:
                list_of_values.append(quantity.as_vector_with_length(1))
            instance._set_values([key], self.attribute_names, list_of_values)
            
    
    
    class CalculatedAttribute:
        
        def  __init__(self, function):
            self.function = function
            arguments, varargs, kwargs, defaults = inspect.getargspec(function)
            self.attribute_names = arguments 
        
        def _get_values(self, instance):
            values = instance._get_values(instance._get_keys(), self.attribute_names)
            return self.function(*values)
        
        def get_value_for_particle(self, instance,  key):
            values = instance._get_values([key], self.attribute_names)
            return self.function(*values)[0]
              

        
            
    def __init__(self):
        object.__setattr__(self, "_vector_attributes", self.GLOBAL_VECTOR_ATTRIBUTES.copy())
        object.__setattr__(self, "_calculated_attributes", self.GLOBAL_CALCULATED_ATTRIBUTES.copy())
        object.__setattr__(self, "_private", self.PrivateProperties())
    
    def __getattr__(self, name_of_the_attribute):
        if name_of_the_attribute == 'key':
            return self._get_keys()
        elif name_of_the_attribute in self._vector_attributes:
            return self._vector_attributes[name_of_the_attribute]._get_values(self)
        elif name_of_the_attribute in self._calculated_attributes:
            return self._calculated_attributes[name_of_the_attribute]._get_values(self)
        else:
            return self._get_values(self._get_keys(), [name_of_the_attribute])[0]
    
    def __setattr__(self, name_of_the_attribute, value):
        if name_of_the_attribute in self._vector_attributes:
            self._vector_attributes[name_of_the_attribute]._set_values(self, value)
        else:
            self._set_values(self._get_keys(), [name_of_the_attribute], [value])
    
    def _get_value_of_attribute(self, key, attribute):
        if attribute in self._vector_attributes:
            return self._vector_attributes[attribute].get_value_for_particle(self, key)
        elif attribute in self._calculated_attributes:
            return self._calculated_attributes[attribute].get_value_for_particle(self, key)
        else:
            return self._get_values([key], [attribute])[0][0]
            
    
    def _set_value_of_attribute(self, key, attribute, value):
        if attribute in self._vector_attributes:
            return self._vector_attributes[attribute].set_value_for_particle(self, key, value)
        else:
            return self._set_values([key], [attribute], value.as_vector_with_length(1))
        
    #
    # Particle storage interface
    #
    
    def _get_values(self, keys, attributes):
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
    
    def _get_attributes(self):
        return []
    
    def _get_state_attributes(self):
        return []
    
    def _real_particles(self):
        return self
    #
    #
    #
    
    def _values_of_particle(self, key):
        attributes = self._get_attributes()
        keys = [key]
        values = self._get_values(keys, attributes)
        
        for attribute, val in zip(attributes, values):
            yield attribute, val[0]
    
    def add_vector_attribute(self, name_of_the_attribute, name_of_the_components):
        """
        Define a vector attribute, coupling two or more scalar attributes into
        one vector attribute. 
        
        :argument name_of_the_attribute: Name to reference the vector attribute by. 
        :argument name_of_the_components: List of strings, each string a name of a scalar attribute.
        
        >>> particles = Particles(2)
        >>> particles.vx = [1.0 , 2.0] | units.m / units.s
        >>> particles.vy = [3.0 , 4.0] | units.m / units.s
        >>> particles.add_vector_attribute("v", ["vx","vy"])
        >>> particles.v
        quantity<[[ 1.  3.], [ 2.  4.]] m / s>
        
        """
        
        self._vector_attributes[name_of_the_attribute] = self.VectorAttribute(name_of_the_components)
    
    @classmethod
    def add_global_vector_attribute(cls, name_of_the_attribute, name_of_the_components):
        """
        Define a *global* vector attribute, coupling two or more scalar attributes into
        one vector attribute. The vector will be defined for all particle sets
        created after calling this function.
        
        :argument name_of_the_attribute: Name to reference the vector attribute by. 
        :argument name_of_the_components: List of strings, each string a name of a scalar attribute.
        
        
        >>> Particles.add_global_vector_attribute('vel', ['vx','vy'])
        >>> particles = Particles(2)
        >>> particles.vx = [1.0 , 2.0] | units.m / units.s
        >>> particles.vy = [3.0 , 4.0] | units.m / units.s
        >>> particles.vel
        quantity<[[ 1.  3.], [ 2.  4.]] m / s>
        
        """
        cls.GLOBAL_VECTOR_ATTRIBUTES[name_of_the_attribute] = cls.VectorAttribute(name_of_the_components)
    
    
    def add_calculated_attribute(self, name_of_the_attribute, function):
        """
        Define a read-only calculated attribute, values for the attribute are
        calculated using the given function. The functions argument 
        names are interperted as attribute names. For example, if
        the given function is::
        
            def norm(x, y):
                return (x*x + y*y).sqrt()
                
        The attributes "x" and "y" will be retrieved from the particles
        and send to the "norm" function.
        
        The calculated values are not stored on the particles. Values
        are recalculated every time this attribute is accessed.
        
        :argument name_of_the_attribute: Name to reference the attribute by. 
        :argument function: Function to call, when attribute is accessed
          
        
        >>> particles = Particles(2)
        >>> particles.x = [1.0 , 2.0] | units.m 
        >>> particles.y = [3.0 , 4.0] | units.m
        >>> particles.add_calculated_attribute("xy", lambda x, y : x * y)
        >>> print particles.xy
        [3.0, 8.0] m**2
        >>> particles[0].x = 4.0 | units.m
        >>> print particles.xy
        [12.0, 8.0] m**2
        >>> print particles[1].xy
        8.0 m**2
        
        """
        
        self._calculated_attributes[name_of_the_attribute] = self.CalculatedAttribute(function)
    
    
    @classmethod
    def add_global_calculated_attribute(cls, name_of_the_attribute, function):
        """
        Define a *global* vector attribute, coupling two or more scalar attributes into
        one vector attribute. The vector will be defined for all particle sets
        created after calling this function.
        
        :argument name_of_the_attribute: Name to reference the vector attribute by. 
        :argument name_of_the_components: List of strings, each string a name of a scalar attribute.
        
        
        >>> Particles.add_global_calculated_attribute("xy", lambda x, y : x * y)
        >>> particles = Particles(2)
        >>> particles.x = [1.0 , 2.0] | units.m 
        >>> particles.y = [3.0 , 4.0] | units.m
        >>> print particles.xy
        [3.0, 8.0] m**2
        
        """
        cls.GLOBAL_CALCULATED_ATTRIBUTES[name_of_the_attribute] = cls.CalculatedAttribute(function)
    
    
    #
    # public API
    #
    def __iter__(self):
        index = 0
        for key in self._get_keys():
            p = Particle(key, self._real_particles())
            yield p
            index += 1

    def __getitem__(self, index):
        return Particle(self._get_keys()[index], self._real_particles())
        
    def __len__(self):
        return len(self._get_keys())
        
     
    def copy(self):
        """
        Creates a new in memory particle set and copies
        all attributes and values into this set. The history
        of the set is not copied over.
        """
        attributes = self._get_attributes()
        keys = self._get_keys()
        values = self._get_values(None, attributes)
        result = Particles()
        result._set_particles(keys, attributes, values)
        result._private._vector_attributes = self._private._vector_attributes.copy()
        return result
         
    
    def copy_values_of_attribute_to(self, attribute_name, particles):
        """
        Copy values of the attributes from this set to the 
        other set. Will only copy values for the particles
        in both sets
        """
        channel = self.new_channel_to(particles)
        channel.copy_attributes([attribute_name])
    
    def new_channel_to(self, other):
        return ParticleInformationChannel(self, other)
        
    def add_particles(self, particles):
        """
        Adds particles from the supplied set to this set. Attributes
        and values are copied over. 
        
        **Note** For performance reasons the particles
        are not checked for duplicates. When the same particle 
        is part of both sets errors may occur.
        
        :parameter particles: set of particles to copy values from
        
        >>> particles1 = Particles(2)
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles2 = Particles(2)
        >>> particles2.x = [3.0, 4.0] | units.m
        >>> particles1.add_particles(particles2)  # doctest:+ELLIPSIS
        <amuse.support.data.core.ParticlesSubset object at 0x...>
        >>> print len(particles1)
        4
        >>> print particles1.x
        [1.0, 2.0, 3.0, 4.0] m
        """
        attributes = particles._get_attributes()
        keys = particles._get_keys()
        values = particles._get_values(keys, attributes)
        self._set_particles(keys, attributes, values)
        return ParticlesSubset(self._real_particles(), keys)
    
    
    def add_particle(self, particle):
        """
        Add one particle to the set. 
        
        :parameter particle: particle to add
        
        >>> particles = Particles()
        >>> print len(particles)
        0
        >>> particle = Particle()
        >>> particle.x = 1.0 | units.m
        >>> particles.add_particle(particle)  # doctest:+ELLIPSIS
        <amuse.support.data.core.Particle object at ...>
        >>> print len(particles)
        1
        >>> print particles.x
        [1.0] m
        
        """
        return self.add_particles(particle.as_set())[0]
        
    
    def remove_particles(self, particles):
        """
        Removes particles from the supplied set from this set.
        
        :parameter particles: set of particles to remove from this set
        
        >>> particles1 = Particles(2)
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles2 = Particles()
        >>> particles2.add_particle(particles1[0]) # doctest:+ELLIPSIS
        <amuse.support.data.core.Particle object at ...>
        >>> particles1.remove_particles(particles2)
        >>> print len(particles1)
        1
        >>> print particles1.x
        [2.0] m
        """
        keys = particles._get_keys()
        self._remove_particles(keys)
        
    
    def remove_particle(self, particle):
        """
        Removes a particle from this set.
        
        Result is undefined if particle is not part of the set
        
        :parameter particle: particle to remove from this set
        
        >>> particles1 = Particles(2)
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles1.remove_particle(particles1[0])
        >>> print len(particles1)
        1
        >>> print particles1.x
        [2.0] m
        """
        self.remove_particles(particle.as_set())
        
    def synchronize_to(self, other_particles):
        """
        Synchronize the particles of this set
        with the contents of the provided set.
        After this call the proveded set will have
        the same particles as the given set. This call
        will check if particles have been removed or
        added it will not copy values of existing particles
        over.
        
        :parameter other_particles: particle set wich has to be updated
        
        >>> particles = Particles(2)
        >>> particles.x = [1.0, 2.0] | units.m
        >>> copy = particles.copy()
        >>> new_particle = Particle()
        >>> new_particle.x = 3.0 | units.m
        >>> particles.add_particle(new_particle)# doctest:+ELLIPSIS
        <amuse.support.data.core.Particle object at ...>
        >>> print particles.x
        [1.0, 2.0, 3.0] m
        >>> print copy.x
        [1.0, 2.0] m
        >>> particles.synchronize_to(copy)
        >>> print copy.x
        [1.0, 2.0, 3.0] m
        
        """
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
        channel.copy_attributes(self._get_state_attributes())   
    
    def to_set(self):
        """
        Returns a subset view on this set. The subset
        will contain all particles of this set.
        """
        return ParticlesSubset(self._real_particles(), self._get_keys())
    
    
    def select(self, selection_function, attributes):
        """
        Returns a subset view on this set. The subset
        will contain all particles for which the selection
        function returned True. The selection function 
        is called with scalar quantities defined by
        the attributes parameter
        
        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> particles.x = [1.0, 2.0, 3.0] | units.m
        >>> subset = particles.select(lambda x : x > 15.0 | units.kg, ["mass"])
        >>> print subset.mass
        [20.0, 30.0] kg
        >>> print subset.x
        [2.0, 3.0] m
        
        """
        keys = self._get_keys()
        values = self._get_values(keys, attributes)
        selected_keys = []
        for index in range(len(keys)):
            key = keys[index]
            arguments = [None] * len(attributes)
            for attr_index, attribute in enumerate(attributes):
                arguments[attr_index] = values[attr_index][index]
            if selection_function(*arguments):
                selected_keys.append(key)
        return ParticlesSubset(self._real_particles(), selected_keys)
    
    def difference(self, other):
        """
        Returns a new subset containing the difference between
        this set and the provided set.
        
        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> particles.x = [1.0, 2.0, 3.0] | units.m
        >>> subset = particles.select(lambda x : x > 15.0 | units.kg, ["mass"])
        >>> less_than_15kg = particles.difference(subset)
        >>> len(subset)
        2
        >>> len(less_than_15kg)
        1
        
        """
        return self.to_set().difference(other)
        
    def has_duplicates(self):
        """
        Returns two when a set contains a particle with the
        same key more than once.
        """
        return len(self) != len(set(self._get_keys()))
    
    
    def as_subset_in(self, other):
        selected_keys = filter(lambda x : other._has_key(x), self._get_keys())
        return other._subset(selected_keys)
        
    def _subset(self, keys):
        return ParticlesSubset(self._real_particles(), keys)
        
        
    def __dir__(self):
        """
        Utility function for introspection of paricle objects
        
        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> particles.x = [1.0, 2.0, 3.0] | units.m
        >>> print 'mass' in dir(particles)
        True
        >>> print 'x' in dir(particles)
        True
        
        """
        result = []
        result.extend(dir(type(self)))
        result.extend(self._attributes_for_dir())
        return result

    
    def _attributes_for_dir(self):
        result = []
        result.extend(self._get_attributes())
        result.extend(self._calculated_attributes.keys())
        result.extend(self._vector_attributes.keys())
        return result
        
class Particles(AbstractParticleSet):
    """
    A set of particles. Attributes and values are stored in
    a private storage model. This storage model can store
    the values in the python memory space, in the memory space
    of the code or in a HDF5 file. By default the storage
    model is in memory.
    
    
    
    """
    def __init__(self, size = 0, storage = None):
        AbstractParticleSet.__init__(self)
        
        if storage is None:
            self._private.attribute_storage = InMemoryAttributeStorage()
        else:
            self._private.attribute_storage = storage
    
        if size > 0:
            particle_keys = UniqueKeyGenerator.next_set_of_keys(size)
            self._set_particles(particle_keys)
            
        self._private.previous = None
        self._private.timestamp = None
        
    def savepoint(self, timestamp=None):
        instance = type(self)()
        instance._private.attribute_storage = self._private.attribute_storage.copy()
        instance._private.timestamp = timestamp
        instance._private.previous = self._private.previous
        self._private.previous = instance
    
    def get_timestamp(self):
        return self._private.timestamp
        
    def iter_history(self):
        current = self
        while not current is None:
            yield current
            current = current._private.previous
    
    def previous_state(self):
        return self._private.previous
        
    @property
    def history(self):
        return reversed(list(self.iter_history()))
        
    def get_timeline_of_attribute(self, particle_key, attribute):
        print attribute
        timeline = []
        for x in self.history:
            if x._has_key(particle_key):
                timeline.append((x._private.timestamp, x._get_value_of_attribute(particle_key, attribute)))
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
        
    def _get_state_attributes(self):
        return self._private.attribute_storage._state_attributes()
        
        
    def _get_keys(self):
        return self._private.attribute_storage._get_keys()
        
    def _has_key(self, key):
        return self._private.attribute_storage._has_key(key)
        
    
    

class ParticlesSubset(AbstractParticleSet):
    """A subset of particles. Attribute values are not
    stored by the subset. The subset provides a limited view
    to the particles. 
    
    Particle subset objects are not supposed to be created
    directly. Instead use the ``to_set`` or ``select`` methods.
    """
    def __init__(self, particles, keys):
        AbstractParticleSet.__init__(self)
        
        self._private.particles = particles
        self._private.keys = numpy.array(keys)
        self._private.set_of_keys = set(keys)
              
        
    def _set_particles(self, keys, attributes = [], values = []):
        """
        Adds particles from to the subset, also
        adds the particles to the superset
        """
        self._private.keys = numpy.concatenate((self.keys,  numpy.array(keys)))
        self._private.set_of_keys += set(keys)
        self._private.particles._set_particles(keys, attributes, values)
        
    def _remove_particles(self, keys):
        """
        Removes particles from the subset, does not remove particles
        from the super set
        """
        set_to_remove = set(keys)
        self._private.set_of_keys -= set_to_remove
        index = 0
        indices = []
        for x in self._private.keys:
            if x in set_to_remove:
                indices.append(index)
            index += 1
        self._private.keys =  numpy.delete(self._private.keys,indices)
    
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
    
    def _get_state_attributes(self):
        return self._private.particles._get_state_attributes(self)
            
        
        
    def _real_particles(self):
        return self._private.particles
        
    def difference(self, other):
        new_set_of_keys = self._private.set_of_keys.difference(other.to_set()._private.set_of_keys)
        return ParticlesSubset(self._private.particles, list(new_set_of_keys))
        
    def union(self, other):
        """
        Returns a new subset containing the union between
        this set and the provided set.
        
        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> subset1 = particles.select(lambda x : x > 25.0 | units.kg, ["mass"])
        >>> subset2 = particles.select(lambda x : x < 15.0 | units.kg, ["mass"])
        >>> union = subset1.union(subset2)
        >>> len(union)
        2
        >>> sorted(union.mass.value_in(units.kg))
        [10.0, 30.0]
        """
        
        new_set_of_keys = self._private.set_of_keys.union(other.to_set()._private.set_of_keys)
        return ParticlesSubset(self._private.particles, list(new_set_of_keys))
    
    def to_set(self):
        return self



class ParticlesWithUnitsConverted(AbstractParticleSet):
    """
    A view on a particle sets. Used when to convert
    values between incompatible sets of units. For example to
    convert from si units to nbody units.
    
    The converter must have implement the ConverterInterface.
    
    >>> from amuse.support.units import nbody_system
    >>> particles_nbody = Particles(2)
    >>> particles_nbody.x = [10.0 , 20.0] | nbody_system.length
    >>> convert_nbody = nbody_system.nbody_to_si(10 | units.kg , 5 | units.m )
    >>> particles_si = ParticlesWithUnitsConverted(
    ...     particles_nbody,
    ...     convert_nbody.as_converter_from_si_to_nbody())
    ...
    >>> print particles_nbody.x
    [10.0, 20.0] nbody length
    >>> print particles_si.x
    [50.0, 100.0] m
    >>> particles_si.x = [200.0, 400.0] | units.m
    >>> print particles_nbody.x
    [40.0, 80.0] nbody length

    
    """
    
    class ConverterInterface(object):
        """
        Interface definition for the converter.
        
        source
            The source quantity is in the units of the user of a
            ParticlesWithUnitsConverted object
        target
            The target quantity must be in the units of the
            internal particles set.
        
        """
        def from_source_to_target(quantity):
            """
            Converts the quantity from the source units
            to the target units.
            
            :parameter quantity: quantity to convert
            """
            return quantity
            
        
        def from_target_to_source(quantity):
            """
            Converts the quantity from the target units
            to the source units.
            
            :parameter quantity: quantity to convert
            """
            return quantity
            
        
    def __init__(self, particles, converter):
        AbstractParticleSet.__init__(self)
        
        self._private.particles = particles
        self._private.converter = converter
              
        
    def _set_particles(self, keys, attributes = [], values = []):
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_source_to_target(quantity)
            converted_values.append(converted_quantity)
        self._private.particles._set_particles(keys, attributes, converted_values)
        
    def _remove_particles(self, keys):
        self._private.particles._remove_particles(keys)
        
    def _get_values(self, keys, attributes):
        values = self._private.particles._get_values(keys, attributes)
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_target_to_source(quantity)
            converted_values.append(converted_quantity)
        return converted_values
        
    def _set_values(self, keys, attributes, values):
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_source_to_target(quantity)
            converted_values.append(converted_quantity)
        self._private.particles._set_values(keys, attributes, converted_values)
    
    def _get_attributes(self):
        return self._private.particles._get_attributes()
        
    def _get_state_attributes(self):
        return self._private.particles._get_state_attributes()
        
    def _get_keys(self):
        return self._private.particles._get_keys()
        
    def _has_key(self, key):
        return self._private.particles._has_key(key)
        
    def to_set(self):
        return ParticlesSubset(self, self._get_keys())
    
            

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
        self._reindex()
        data = self.from_particles._get_values(self.keys, attributes)
        self.to_particles._set_values(self.keys, attributes, data)
    
    def copy(self):
        self._reindex()
        self.copy_attributes(self.from_particles._get_attributes())
    

        
class Stars(Particles):

    def __init__(self, size = 0):
        Particles.__init__(self, size)
    
    
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

class Particle(object):
    """A physical object or a physical region simulated as a 
    physical object (cloud particle).
    
    All attributes defined on a particle are specific for 
    that particle (for example mass or position). A particle contains 
    a set of attributes, some attributes are *generic* and applicable
    for multiple modules. Other attributes are *specific* and are 
    only applicable for a single module.
    
    
    """
    __slots__ = ["key", "particles_set"]
    
    
    
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
       
        if isinstance(new_value_for_the_attribute, values.Quantity):
            self.particles_set._set_value_of_attribute(self.key, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise Exception("attribute "+name_of_the_attribute+" does not have a valid value, values must have a unit")
    
    def __getattr__(self, name_of_the_attribute):
         return self.particles_set._get_value_of_attribute(self.key, name_of_the_attribute)
    
    def children(self):
        return self.particles_set.select(lambda x : x == self.key | units.none, ["parent_id"])
    
    def descendents(self):
        result = self.children()
        stack = list(result)
        while len(stack) > 0:
            current = stack.pop()
            children = current.children()
            result = result.union(children)
            stack.extend(children)
        return result
        
    def add_child(self, child):
        if self.particles_set != child.particles_set:
            raise Exception("The parent and child particles should be in the same set")
        
        child.parent_id = self.key | units.none
                
    def __str__(self):
        output = 'Particle '
        output += str(self.key)
        output += ''
        output += '\n'
        for name, value in self.particles_set._values_of_particle(self.key):
            output += name
            output += ': {'
            output += str(value)
            output += '}, '
            output += '\n'
        return output
    
    def __dir__(self):
        result = []
        result.extend(dir(type(self)))
        result.extend(self.particles_set._attributes_for_dir())
        return result
        
        
    def set_default(self, attribute, quantity):
        if not attribute in self.particles_set._get_attributes():
            self.particles_set._set_value_of_attribute(self, attribute, quantity)
            
    def get_timeline_of_attribute(self, attribute):
        return self.particles_set.get_timeline_of_attribute(self.key, attribute)
        
    def as_set(self):
        """
        Returns a subset view on the set containing this particle. The
        subset view includes this particle and no other particles.
        
        >>> particles = Particles(2)
        >>> particles.x = [1.0, 2.0] | units.m
        >>> particle2 = particles[1]
        >>> print particle2.x
        2.0 m
        >>> particles_with_one_particle = particle2.as_set()
        >>> len(particles_with_one_particle)
        1
        >>> print particles_with_one_particle.x
        [2.0] m
        """
        return ParticlesSubset(self.particles_set, [self.key])
        

AbstractParticleSet.add_global_vector_attribute("position", ["x","y","z"])
AbstractParticleSet.add_global_vector_attribute("velocity", ["vx","vy","vz"])


class InterfaceWithUnitsConverted(object):
    class UnitsConvertionMethod(object):
        def __init__(self, real_method, converter):
            self.real_method = real_method
            self.converter = converter
            
        def __call__(self, *list_arguments, **keyword_arguments):
            converted_list_arguments = [self.from_source_to_target(x) for x in list_arguments]
            converted_keyword_arguments = {}
            for key, value in keyword_arguments:
                converted_keyword_arguments[key] = self.from_source_to_target(value)
            
            result = self.real_method(*converted_list_arguments, **converted_keyword_arguments)
            return self.from_target_to_source(result)
            
        def from_source_to_target(self, x):
            if isinstance(x, values.Quantity):
                return self.converter.from_source_to_target(quantity)
            else:
                x
                
        def from_target_to_source(self, x):
            if isinstance(x, values.Quantity):
                return self.converter.from_target_to_source(quantity)
            else:
                x
    
    def __init__(self, interface, converter):
        self.converter = converter
        self.interface = interface
        
    def __getattr__(self, name):
        attribute = getattr(self.interace, name)
        if inspect.ismethod(attribute):
            return self.UnitsConvertionMethod(attribute)
        elif isinstance(attribute, AbstractParticleSet):
            return ParticlesWithUnitsConverted(attribute, self.converter)
        elif isinstance(attribute, values.Quantity):
            return self.converter.from_target_to_source(quantity)
        else:
            return attribute
            
    def __dir__(self):
        return dir(self.interface)
        
        
