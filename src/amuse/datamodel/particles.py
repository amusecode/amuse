from amuse.support.core import CompositeDictionary
from amuse.support import exceptions
from amuse.datamodel.base import *
from amuse.datamodel.memory_storage import *
from amuse.units import constants
from amuse.units import units
from amuse.units import quantities
from amuse.units.quantities import Quantity
from amuse.units.quantities import new_quantity
from amuse.units.quantities import zero
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.units.quantities import AdaptingVectorQuantity

import numpy

class AbstractParticleSet(AbstractSet):
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
    
    
    Particle sets can be sliced.
    
    >>> particles = Particles(3)
    >>> particles.x = [10.0, 20.0, 30.0] | units.m
    >>> particles[1:].x
    quantity<[20.0, 30.0] m>
    
    
    Particle sets can be copied.
    
    >>> particles = Particles(3)
    >>> particles.x = [10.0, 20.0, 30.0] | units.m
    >>> copy = particles.copy()
    >>> particles.x = 2.0 | units.m
    >>> particles.x
    quantity<[2.0, 2.0, 2.0] m>
    >>> copy.x
    quantity<[10.0, 20.0, 30.0] m>
    
    
    Particle sets can be added together. 
    Attribute values are not stored by the resulting subset. The subset
    provides a view on two or more sets of particles. Changing attributes
    of the sum of sets will also change the attributes of each original
    set, contrary to copy().
    
    >>> particles = Particles(4)
    >>> particles1 = particles[:2]
    >>> particles1.x = [1.0, 2.0] | units.m
    >>> particles2 = particles[2:]
    >>> particles2.x = [3.0, 4.0] | units.m
    >>> new_set = particles1 + particles2
    >>> print len(new_set)
    4
    >>> print new_set.x
    [1.0, 2.0, 3.0, 4.0] m
    
    
    Particle sets can be subtracted from each other. 
    Like with addition, attribute values are not stored by the resulting
    subset.
    
    >>> particles = Particles(4)
    >>> particles.x = [1.0, 2.0, 3.0, 4.0] | units.m
    >>> junk = particles[2:]
    >>> new_set = particles - junk
    >>> print len(new_set)
    2
    >>> print new_set.x
    [1.0, 2.0] m
    >>> print particles.x
    [1.0, 2.0, 3.0, 4.0] m
    
    
    Particle sets can have instance based or global vector attributes.
    A particle set stores a list of scalar values for each attribute.
    Some attributes are more naturally accessed as lists 
    of vector values. Once defined, a particle set can
    convert the scalar values of 2 or more attributes into one
    vector attribute.
    
    >>> from amuse.support.data import particle_attributes
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
    
    
    GLOBAL_DERIVED_ATTRIBUTES = {}
    
    
    def __init__(self, original = None):
        AbstractSet.__init__(self, original)
    
    
    def check_attribute(self, value):
        if not (isinstance(value, Quantity) or isinstance(value, Particle) or isinstance(value, AbstractParticleSet)):
            if hasattr(value, "__iter__"):
                result = AdaptingVectorQuantity()
                for subvalue in value:
                    result.append(self.check_attribute(subvalue))
                return result
            raise AttributeError("Can only assign quantities or other particles to an attribute.")
        return value
            
    #
    # Particle storage interface
    #
    
    def add_particles_to_store(self, keys, attributes, values):
        pass
        
    def remove_particles_from_store(self, keys):
        pass
        
    def get_values_in_store(self, keys, attributes):
        pass
    
    def set_values_in_store(self, keys, attributes, values):
        pass
        
    def get_all_keys_in_store(self):
        return []
        
    def has_key_in_store(self):
        return False
    
    def get_attribute_names_defined_in_store(self):
        return []
    
        
    #
    #
    #
    
    def _values_of_particle(self, key):
        attributes = self.get_attribute_names_defined_in_store()
        keys = [key]
        values = self.get_values_in_store(keys, attributes)
        
        for attribute, val in zip(attributes, values):
            yield attribute, val[0]
        
    #
    # public API
    #
    def __iter__(self):
        original_set = self._original_set()
        for key in self.get_all_keys_in_store():
            yield Particle(key,original_set)

    def __getitem__(self, index):
        keys = self.get_all_keys_in_store()[index]
        
        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Particle(self.get_all_keys_in_store()[index], self._original_set())
    
    def get_all_particles_at(self, *indices):
        all_keys = self.get_all_keys_in_store()
        selectedkeyes = [all_keys[x] for x in indices]
        return self._subset(selectedkeyes)
        
    def __str__(self):
        """
        Display string of a particle set.
        
        >>> p0 = Particle(10)
        >>> p1 = Particle(11)
        >>> particles = Particles()
        >>> particles.add_particle(p0) # doctest:+ELLIPSIS
        <amuse.datamodel.particles.Particle object at ...>
        >>> particles.add_particle(p1) # doctest:+ELLIPSIS
        <amuse.datamodel.particles.Particle object at ...>
        >>> particles.x = [4.0 , 3.0] | units.m
        >>> particles.y = [5.0 , 2.0] | units.km
        >>> print particles 
                         key            x            y
                           -            m           km
        ====================  ===========  ===========
                          10    4.000e+00    5.000e+00
                          11    3.000e+00    2.000e+00
        ====================  ===========  ===========

        """
        return self.to_string()
        
        
    def to_string(self, attributes_to_show = None):
        """
        Display string of a particle set.
        
        >>> p0 = Particle(10)
        >>> p1 = Particle(11)
        >>> particles = Particles()
        >>> particles.add_particle(p0) # doctest:+ELLIPSIS
        <amuse.datamodel.particles.Particle object at ...>
        >>> particles.add_particle(p1) # doctest:+ELLIPSIS
        <amuse.datamodel.particles.Particle object at ...>
        >>> particles.x = [4.0 , 3.0] | units.m
        >>> particles.y = [5.0 , 2.0] | units.km
        >>> print particles 
                         key            x            y
                           -            m           km
        ====================  ===========  ===========
                          10    4.000e+00    5.000e+00
                          11    3.000e+00    2.000e+00
        ====================  ===========  ===========

        """
        attributes = sorted(self.get_attribute_names_defined_in_store())
        if attributes_to_show:
            attributes = [x for x in attributes if x in attributes_to_show]
        
        format_float = '{0: >11.3e}'.format
        format_str20 = '{0: >20}'.format
        format_str11 = '{0: >11}'.format

        columns = map(lambda x : [format_str11(x)], attributes)
        columns.insert(0,[format_str20('key')])
        
        all_values = self.get_values_in_store(self.get_all_keys_in_store(), attributes)
        for index, quantity in enumerate(all_values):
            column = columns[index + 1]
            if hasattr(quantity, 'unit'):
                column.append(format_str11(str(quantity.unit)))
            else:
                column.append('none')
                
            column.append('=' * 11)
            if len(quantity) > 40:
                if hasattr(quantity, 'unit'):
                    if numpy.issubdtype(quantity.number.dtype, float):
                        values_to_show = list(map(format_float,quantity.number[:20]))
                        values_to_show.append(format_str11('...'))
                        values_to_show.extend(map(format_float,quantity.number[-20:]))
                    else:
                        values_to_show = list(map(format_str11,quantity.number[:20]))
                        values_to_show.append(format_str11('...'))
                        values_to_show.extend(map(format_str11,quantity.number[-20:]))
                else:
                    values_to_show = list(map(format_str11,quantity[:20]))
                    values_to_show.append(format_str11('...'))
                    values_to_show.extend(map(format_str11,quantity[-20:]))
            else:
                if hasattr(quantity, 'unit'):
                    if numpy.issubdtype(quantity.number.dtype, float):
                        values_to_show = map(format_float,quantity.number)
                    else:
                        values_to_show = map(format_str11,quantity.number)
                else:
                    values_to_show = map(format_str11,quantity)
                    
            column.extend(values_to_show)
            column.append('=' * 11)
            
        column = columns[0]
        column.append(format_str20("-"))
        column.append('=' * 20)
        particle_keys = self.get_all_keys_in_store()
        if len(particle_keys) > 40:
            values_to_show = list(map(format_str20, particle_keys[:20]))
            values_to_show.append(format_str20('...'))
            values_to_show.extend(map(format_str20, particle_keys[-20:]))
        else:
            values_to_show = map(format_str20,particle_keys)
                    
        column.extend(values_to_show)
            
        column.append('=' * 20)
        
        rows = []
        for i in range(len(columns[0])):
        
            row = [x[i] for x in columns]
            rows.append(row)
            
        lines = map(lambda  x : '  '.join(x), rows)
        return '\n'.join(lines)    
        
    def _get_particle(self, key):
        if self.has_key_in_store(key):
            return Particle(key, self._original_set())
        else:
            return None
            
    def can_extend_attributes(self):
        return self._original_set().can_extend_attributes()
        
        
    def copy(self):
        """
        Creates a new in memory particle set and copies
        all attributes and values into this set. The history
        of the set is not copied over.
        """
        attributes = self.get_attribute_names_defined_in_store()
        keys = self.get_all_keys_in_store()
        values = self.get_values_in_store(keys, attributes)
        result = self._factory_for_new_collection()()
        result.add_particles_to_store(keys, attributes, values)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
       
        return result
    
    def _is_superset(self):
        return False

    
    def copy_to_memory(self):
        attributes = self.get_attribute_names_defined_in_store()
        keys = self.get_all_keys_in_store()
        values = self.get_values_in_store(keys, attributes)
        result = Particles()
        result.add_particles_to_store(keys, attributes, values)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
       
        return result
    
    def copy_values_of_attribute_to(self, attribute_name, particles):
        """
        Copy values of one attribute from this set to the 
        other set. Will only copy values for the particles
        in both sets. See also :meth:`synchronize_to`.
        
        If you need to do this a lot, setup a dedicated
        channel.
        
        >>> particles1 = Particles(2)
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles2 = particles1.copy()
        >>> print particles2.x
        [1.0, 2.0] m
        >>> p3 = particles1.add_particle(Particle())
        >>> particles1.x = [3.0, 4.0, 5.0] | units.m
        >>> particles1.copy_values_of_attribute_to("x", particles2)
        >>> print particles2.x
        [3.0, 4.0] m
        
        """
        channel = self.new_channel_to(particles)
        channel.copy_attributes([attribute_name])
    
    def new_channel_to(self, other):
        return ParticleInformationChannel(self, other)
        
    def __add__(self, particles):
        """
        Returns a particle subset, composed of the given
        particle(s) and this particle set. Attribute values are
        not stored by the subset. The subset provides a view
        on two or more sets of particles.
        
        :parameter particles: (set of) particle(s) to be added to self.
        
        >>> particles = Particles(4)
        >>> particles1 = particles[:2]
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles2 = particles[2:]
        >>> particles2.x = [3.0, 4.0] | units.m
        >>> new_set = particles1 + particles2
        >>> new_set  # doctest:+ELLIPSIS
        <amuse.datamodel.particles.ParticlesSubset object at 0x...>
        >>> print len(new_set)
        4
        >>> print new_set.x
        [1.0, 2.0, 3.0, 4.0] m
        """
        if isinstance(particles, Particle):
            particles = particles.as_set()
        original_particles_set = self._original_set()
        if set(original_particles_set.key)!=set(particles._original_set().key):
            raise exceptions.AmuseException("Can't create new subset from particles belonging to "
                "separate particle sets. Try creating a superset instead.")
        keys = list(self.key) + list(particles.key)
        new_set = ParticlesSubset(original_particles_set, keys)
        if new_set.has_duplicates():
            raise exceptions.AmuseException("Unable to add a particle, because it was already part of this set.")
        return new_set
    
    def __or__(self, particles):
        """
        Returns a particle superset, composed of the given
        particle(s) and this particle set.
        
        :parameter particles: (set of) particle(s) to be added to self.
         
        >>> particles1 = Particles(2)
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles2 = Particles(2)
        >>> particles2.x = [3.0, 4.0] | units.m
        >>> new_set = particles1 | particles2
        >>> new_set  # doctest:+ELLIPSIS
        <amuse.datamodel.particles.ParticlesSuperset object at 0x...>
        >>> print len(new_set)
        4
        >>> print new_set.x
        [1.0, 2.0, 3.0, 4.0] m
        """
        if isinstance(particles, Particle):
            particles = particles.as_set()
            
        original_particles_set1 = self._original_set()
        original_particles_set2 = particles._original_set()
        
        return ParticlesSuperset((original_particles_set1, original_particles_set2))
    
    def __sub__(self, particles):
        """
        Returns a subset of the set without the given particle(s)
        Attribute values are not stored by the subset. The subset 
        provides a view on two or more sets of particles.
        
        :parameter particles: (set of) particle(s) to be subtracted from self.
        
        >>> particles = Particles(4)
        >>> particles.x = [1.0, 2.0, 3.0, 4.0] | units.m
        >>> junk = particles[2:]
        >>> new_set = particles - junk
        >>> new_set  # doctest:+ELLIPSIS
        <amuse.datamodel.particles.ParticlesSubset object at 0x...>
        >>> print len(new_set)
        2
        >>> print new_set.x
        [1.0, 2.0] m
        >>> print particles.x
        [1.0, 2.0, 3.0, 4.0] m
        """
        if isinstance(particles, Particle):
            particles = particles.as_set()
        new_keys = []
        new_keys.extend(self.get_all_keys_in_store())
        subtract_keys = particles.get_all_keys_in_store()
        for key in subtract_keys:
            if key in new_keys:
                new_keys.remove(key)
            else:
                raise exceptions.AmuseException("Unable to subtract a particle, because "
                    "it is not part of this set.")
        return self._subset(new_keys)
    
    def add_particles(self, particles):
        """
        Adds particles from the supplied set to this set. Attributes
        and values are copied over. 
        
        .. note::
            For performance reasons the particles
            are not checked for duplicates. When the same particle 
            is part of both sets errors may occur.
        
        :parameter particles: set of particles to copy values from
        
        >>> particles1 = Particles(2)
        >>> particles1.x = [1.0, 2.0] | units.m
        >>> particles2 = Particles(2)
        >>> particles2.x = [3.0, 4.0] | units.m
        >>> particles1.add_particles(particles2)  # doctest:+ELLIPSIS
        <amuse.datamodel.particles.ParticlesSubset object at 0x...>
        >>> print len(particles1)
        4
        >>> print particles1.x
        [1.0, 2.0, 3.0, 4.0] m
        """
        attributes = particles.get_attribute_names_defined_in_store()
        keys = particles.get_all_keys_in_store()
        values = particles.get_values_in_store(keys, attributes)
        values = map(self._convert_from_entities_or_quantities, values)
        self.add_particles_to_store(keys, attributes, values)
        return ParticlesSubset(self._original_set(), keys)
    
    
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
        <amuse.datamodel.particles.Particle object at ...>
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
        <amuse.datamodel.particles.Particle object at ...>
        >>> particles1.remove_particles(particles2)
        >>> print len(particles1)
        1
        >>> print particles1.x
        [2.0] m
        """
        keys = particles.get_all_keys_in_store()
        self.remove_particles_from_store(keys)
        
    
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
        
        After this call the `other_particles` set will have
        the same particles as this set.
        
        This call will check if particles have been removed or
        added it will not copy values of existing particles
        over.
        
        :parameter other_particles: particle set wich has to be updated
        
        >>> particles = Particles(2)
        >>> particles.x = [1.0, 2.0] | units.m
        >>> copy = particles.copy()
        >>> new_particle = Particle()
        >>> new_particle.x = 3.0 | units.m
        >>> particles.add_particle(new_particle)# doctest:+ELLIPSIS
        <amuse.datamodel.particles.Particle object at ...>
        >>> print particles.x
        [1.0, 2.0, 3.0] m
        >>> print copy.x
        [1.0, 2.0] m
        >>> particles.synchronize_to(copy)
        >>> print copy.x
        [1.0, 2.0, 3.0] m
        
        """
        other_keys = set(other_particles.get_all_keys_in_store())
        my_keys = set(self.get_all_keys_in_store())
        added_keys = my_keys - other_keys
        removed_keys = other_keys - my_keys
        added_keys = list(added_keys)
        if added_keys:
            attributes = self.get_attribute_names_defined_in_store()
            values = self.get_values_in_store(added_keys, attributes)
            other_particles.add_particles_to_store(added_keys, attributes, values)
        
        removed_keys = list(removed_keys)
        if removed_keys:
            other_particles.remove_particles_from_store(removed_keys)
        
    def copy_values_of_all_attributes_to(self, particles):
        channel = self.new_channel_to(particles)
        channel.copy_attributes(self.get_attribute_names_defined_in_store())   
    
    def as_set(self):
        """
        Returns a subset view on this set. The subset
        will contain all particles of this set.
        
        >>> particles = Particles(3)
        >>> particles.x = [1.0, 2.0, 3.0] | units.m
        >>> subset = particles.as_set()
        >>> print subset.x
        [1.0, 2.0, 3.0] m
        >>> print particles.x
        [1.0, 2.0, 3.0] m
        """
        return self._subset(self.get_all_keys_in_store())
    
    
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
        keys = self.get_all_keys_in_store()
        #values = self._get_values(keys, attributes) #fast but no vectors
        values = map(lambda x: getattr(self, x), attributes)
        selected_keys = []
        for index in range(len(keys)):
            key = keys[index]
            arguments = [None] * len(attributes)
            for attr_index, attribute in enumerate(attributes):
                arguments[attr_index] = values[attr_index][index]
            if selection_function(*arguments):
                selected_keys.append(key)
        return self._subset(selected_keys)
        
    def select_array(self, selection_function, attributes = ()):
        """
        Returns a subset view on this set. The subset
        will contain all particles for which the selection
        function returned True. The selection function 
        is called with a vector quantities containing all
        the values for the attributes parameter.
        
        This function can be faster than the select function
        as it works on entire arrays. The selection_function
        is called once.
        
        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> particles.x = [1.0, 2.0, 3.0] | units.m
        >>> subset = particles.select_array(lambda x : x > 15.0 | units.kg, ["mass"])
        >>> print subset.mass
        [20.0, 30.0] kg
        >>> print subset.x
        [2.0, 3.0] m
        
        
        >>> particles = Particles(1000)
        >>> particles.x = units.m.new_quantity(numpy.arange(1,1000))
        >>> subset = particles.select_array(lambda x : x > (500 | units.m), ("x",) )
        >>> print len(subset)
        499
        """
        keys = self.get_all_keys_in_store()
        #values = self._get_values(keys, attributes) #fast but no vectors
        quantities = [getattr(self, x) for x in attributes]
        
        selections = selection_function(*quantities)
        selected_keys =  numpy.compress(selections, keys)
        
        return self._subset(selected_keys)
    
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
        return self.as_set().difference(other)
        
    def get_timestamp(self):
        return None
    
    def has_duplicates(self):
        """
        Returns True when a set contains a particle with the
        same key more than once. Particles with the same
        key are interpreted as the same particles.
        
        >>> particles = Particles()
        >>> p1 = particles.add_particle(Particle(1))
        >>> p2 = particles.add_particle(Particle(2))
        >>> particles.has_duplicates()
        False
        >>> p3 = particles.add_particle(Particle(1))
        >>> particles.has_duplicates()
        True
        >>> p3 == p1
        True
        """
        return len(self) != len(set(self.get_all_keys_in_store()))
    
    

        
    def _subset(self, keys):
        return ParticlesSubset(self._original_set(), keys)
    
    def reversed(self):
        """
        Returns a subset with the same particles, but with reversed
        sequenctial order (the first particle will become last)
        
        >>> particles = Particles(3)
        >>> particles.radius = [1.0, 2.0, 3.0] | units.m
        >>> r = particles.reversed()
        >>> print r.radius
        [3.0, 2.0, 1.0] m
        
        """
        
        keys = self.get_all_keys_in_store()
        return self._subset(keys[::-1])
        
    def sorted_by_attribute(self, attribute, kind='mergesort'):
        """
        Returns a subset with the same particles, but sorted
        using the given attribute name
        
        :argument: kind, the sort method for supported kinds see
            the numpy.sort documentation
            
        
        >>> particles = Particles(3)
        >>> particles.mass = [2.0, 3.0, 1.0] | units.kg
        >>> particles.radius = [1.0, 2.0, 3.0] | units.m
        >>> sorted = particles.sorted_by_attribute('mass')
        >>> print sorted.mass
        [1.0, 2.0, 3.0] kg
        >>> print sorted.radius
        [3.0, 1.0, 2.0] m
        """
        values = getattr(self, attribute)
        sorted_indices =  values.argsort(kind=kind)
        keys = self.get_all_keys_in_store()
        return self._subset(keys[sorted_indices])
        
        
    def sorted_by_attributes(self, *attributes):
        """
        Returns a subset with the same particles, but sorted
        using the given attribute names. The last attribute name
        in the call is used for the primary sort order, the 
        second-to-last attribute name for the secondary sort order,
        and so on. See also numpy.lexsort
                  
        
        >>> particles = Particles(4)
        >>> particles.mass = [2.0, 3.0, 1.0, 4.0] | units.kg
        >>> particles.radius = [3.0, 2.0, 1.0, 2.0] | units.m
        >>> sorted = particles.sorted_by_attributes('mass', 'radius')
        >>> print sorted.radius
        [1.0, 2.0, 2.0, 3.0] m
        >>> print sorted.mass
        [1.0, 3.0, 4.0, 2.0] kg
        """
        keys = self.get_all_keys_in_store()
        values = self.get_values_in_store(keys, attributes)
        values = [x.number for x in values]
        sorted_indices =  numpy.lexsort(values)
        
        return self._subset(keys[sorted_indices])
    
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
        result.extend(self.get_attribute_names_defined_in_store())
        result.extend(self._derived_attributes.keys())
        return result
        
    
    def all_attributes(self):
        result = []
        result.append('key')
        result.extend(self._attributes_for_dir())
        return result
        

    def is_empty(self):
        return self.__len__()==0
        
    def get_intersecting_subset_in(self, other):
        selected_keys = filter(lambda x : other.has_key_in_store(x), self.get_all_keys_in_store())
        return other._subset(selected_keys)
        
class Particles(AbstractParticleSet):
    """
    A set of particles. Attributes and values are stored in
    a private storage model. This storage model can store
    the values in the python memory space, in the memory space
    of the code or in a HDF5 file. By default the storage
    model is in memory.
    
    
    
    """
    def __init__(self, size = 0, storage = None, keys = None, keys_generator = None, particles = None):
        AbstractParticleSet.__init__(self)
        
        if storage is None:
            self._private.attribute_storage = get_in_memory_attribute_storage_factory()()
        else:
            self._private.attribute_storage = storage
    
        if not particles is None:
            if isinstance(particles,AbstractParticleSet):
                self.add_particles(particles)
            else:
                for x in iter(particles):
                    self.add_particle(x)
        elif size > 0 or not keys is None:
            if keys is None:
                if keys_generator is None:
                    keys_generator = UniqueKeyGenerator
                    
                particle_keys = keys_generator.next_set_of_keys(size)
            else:
                particle_keys = keys
            self.add_particles_to_store(particle_keys)
            
        self._private.previous = None
        self._private.timestamp = None
        
    def savepoint(self, timestamp=None):
        instance = type(self)()
        instance._private.attribute_storage = self._private.attribute_storage.copy()
        instance._private.timestamp = timestamp
        instance._private.previous = self._private.previous
        self._private.previous = instance
        return instance
    
    def get_timestamp(self):
        return self._private.timestamp
        
    def iter_history(self):
        current = self._private.previous
        while not current is None:
            yield current
            current = current._private.previous
            
    
    def get_state_at_timestamp(self, timestamp):
        previous_timestamp = None
        states_and_distances = []
        for state in self.iter_history():
            timestamp_of_state = state.get_timestamp()
            if timestamp_of_state is None:
                continue
            distance = abs(timestamp_of_state - timestamp)
            states_and_distances.append((state, distance,))
        
        if len(states_and_distances) == 0:
            raise exceptions.AmuseException("You asked for a state at timestamp '{0}', but the set does not have any saved states so this state cannot be returned")
            
        accompanying_state, min_distance = states_and_distances[0]
        for state, distance  in states_and_distances:
            if distance < min_distance:
                min_distance = distance
                accompanying_state = state
        
        return accompanying_state
            
    def previous_state(self):
        return self._private.previous
        
    @property
    def history(self):
        return reversed(list(self.iter_history()))
        
    def get_timeline_of_attribute(self, particle_key, attribute):
        timeline = []
        for x in self.history:
            if x.has_key_in_store(particle_key):
                timeline.append((x._private.timestamp, x._get_value_of_attribute(particle_key, attribute)))
        return timeline

    def get_timeline_of_attribute_as_vector(self, particle_key, attribute):
        timeline = AdaptingVectorQuantity()
        chrono_values = AdaptingVectorQuantity()
        for x in self.history:
            if x.has_key_in_store(particle_key):
                timeline.append(x._private.timestamp)
                chrono_values.append(x._get_value_of_attribute(particle_key, attribute))
        return timeline, chrono_values
    
    def get_timeline_of_attributes(self, particle_key, attributes):
        result = map(lambda x: [], range(len(attributes)+1))
        units = map(lambda x: None, range(len(attributes)+1))
        
        for x in self.history:
            if x.has_key_in_store(particle_key):
                if  units[0] is None:
                    units[0] = x._private.timestamp.unit
                for i, attribute in enumerate(attributes):
                    quantity = x._get_value_of_attribute(particle_key, attribute)
                    if  units[i+1] is None:
                        units[i+1] = quantity.unit
                    result[i+1].append(quantity.value_in(units[i+1]))
                    
        return list(map(lambda value,unit : unit.new_quantity(value), result, units))

                    
            
    def add_particles_to_store(self, keys, attributes = [], values = []):
        self._private.attribute_storage.add_particles_to_store(keys, attributes, values)
        
    def remove_particles_from_store(self, keys):
        self._private.attribute_storage.remove_particles_from_store(keys)
    
    def get_values_in_store(self, keys, attributes, indices = None, version = None):
        return self._private.attribute_storage.get_values_in_store(keys, attributes)
        
    def set_values_in_store(self, keys, attributes, values):
        self._private.attribute_storage.set_values_in_store(keys, attributes, values)
    
    def get_attribute_names_defined_in_store(self):
        return self._private.attribute_storage.get_defined_attribute_names()
        
    def get_all_keys_in_store(self):
        return self._private.attribute_storage.get_all_keys_in_store()
        
    def has_key_in_store(self, key):
        return self._private.attribute_storage.has_key_in_store(key)

    def get_value_in_store(self, key, attribute):
        return self._private.attribute_storage.get_value_in_store(key, attribute)

    def _factory_for_new_collection(self):
        return Particles
        
    def can_extend_attributes(self):
        return self._private.attribute_storage.can_extend_attributes()
    
class BoundSupersetParticlesFunctionAttribute(object):
    def  __init__(self, name, superset):
        self.name = name
        self.superset = superset
        self.subsetfunctions = []
        
    def add_subsetfunction(self, callable):
        self.subsetfunctions.append(callable)
        
    def __call__(self, *list_arguments, **keyword_arguments):
        result = None
        offset = 0
        for x in self.subsetfunctions:
            subset_result = x(*list_arguments, **keyword_arguments)
            if subset_result is None and result is None:
                result = None
            elif hasattr(subset_result, 'unit'):
                if result is None:
                    result = VectorQuantity.zeros(len(self.superset), subset_result.unit)
                    offset = 0
                result[offset:len(subset_result)+offset] = subset_result
                offset += len(subset_result)
            else:
                if result is None:
                    result = []
                result.extend(subset_result)
        return result

class DerivedSupersetAttribute(DerivedAttribute):
                        
            
    def __init__(self, name):
        self.name = name
    
    def get_values_for_entities(self, superset):
        result = None
        offset = 0
        for subset in superset._private.particle_sets:
            
            subset_result =  getattr(subset, self.name)
            if hasattr(subset_result, '__call__'):
                if len(subset) > 0:
                    if result is None:
                        result = BoundSupersetParticlesFunctionAttribute(self.name, superset)
                    result.add_subsetfunction(subset_result)
                
            elif hasattr(subset_result, 'unit'):
                if len(subset_result) == 0:
                    continue
                if result is None:
                    shape = [len(superset),] + list(subset_result.shape[1:])
                    result = VectorQuantity.zeros(shape, subset_result.unit)
                    offset = 0
                try:
                    result[offset:len(subset_result)+offset] = subset_result
                except ValueError:
                    raise AttributeError("Subsets return incompatible quantities for attribute '{0}', attribute cannot be queried from the superset".format(self.name)) 
                offset += len(subset_result)
            else:
                raise exceptions.AmuseException("cannot handle this type of attribute on supersets yet") 
        return result
                
    def set_values_for_entities(self, superset, value):
        raise exceptions.AmuseException("cannot set value of attribute '{0}'")

    def get_value_for_entity(self, superset, key):
        #subset = superset._get_subset_for_key(key)
        #return getattr(subset, self.name)
        raise exceptions.AmuseException("cannot get value of entity for superset")

    def set_value_for_entity(self, superset, key, value):
        #subset = superset._get_subset_for_key(key)
        #return setattr(subset, self.name)
        raise exceptions.AmuseException("cannot get value of entity for superset")

class ParticlesSuperset(AbstractParticleSet):
    """A superset of particles. Attribute values are not
    stored by the superset. The superset provides a view
    on two or more sets of particles. 
    
    Superset objects are not supposed to be created
    directly. Instead use the ``union`` methods.
    
    >>> p1 = Particles(3)
    >>> p1.mass = [10.0, 20.0, 30.0] | units.kg
    >>> p2 = Particles(3)
    >>> p2.mass = [40.0, 50.0, 60.0] | units.kg
    >>> p = ParticlesSuperset([p1, p2])
    >>> print len(p)
    6
    >>> print p.mass
    [10.0, 20.0, 30.0, 40.0, 50.0, 60.0] kg
    >>> p[4].mass = 70 | units.kg
    >>> print p.mass
    [10.0, 20.0, 30.0, 40.0, 70.0, 60.0] kg
    >>> p2[1].mass
    quantity<70.0 kg>
    >>> cp = p.copy()
    >>> print len(cp)
    6
    >>> print cp.mass
    [10.0, 20.0, 30.0, 40.0, 70.0, 60.0] kg
    """
    
    def __init__(self, particle_sets, index_to_default_set=None):
        AbstractParticleSet.__init__(self)
                
        self._private.particle_sets = particle_sets
        self._private.index_to_default_set = index_to_default_set
        
        keys_in_all_subsets = None
        for subset in particle_sets:
            derivedattribute_keys = set(subset._derived_attributes.keys())
            if keys_in_all_subsets is None:
                keys_in_all_subsets = set(derivedattribute_keys)
            else:
                keys_in_all_subsets &= derivedattribute_keys
        
        keys_in_all_subsets -= set(self.GLOBAL_DERIVED_ATTRIBUTES.keys())
        for key in keys_in_all_subsets:
            self._derived_attributes[key] = DerivedSupersetAttribute(key)
        
        if self.has_duplicates():
            raise exceptions.AmuseException("Unable to add a particle, because it was already part of this set.")
        
    
    def can_extend_attributes(self):
        for x in self._private.particle_sets():
            if not x.can_extend_attributes():
                return False
        return True
        
    def __len__(self):
        result = 0
        for set in self._private.particle_sets:
            result += len(set)
            
        return result
    
    def __iter__(self):
        for set in self._private.particle_sets:
            for particle in set:
                yield particle
                
    def _factory_for_new_collection(self):
        return Particles
        
    def __getitem__(self, index):
        offset = 0
        is_index_a_slice = hasattr(index, 'step')
        
        if is_index_a_slice:
            
            for set in self._private.particle_sets:
                length = len(set)
                if index.start < (offset+length):
                    if index.stop > (offset+length):
                        raise Exception('cannot slice accross multiple subsets in a superset')
                    return set[index.start - offset: index.stop - offset: index.step]
                offset += length
            raise Exception('error in interpreting slice on superset')
        else:
            for set in self._private.particle_sets:
                length = len(set)
                if index < (offset+length):
                    return set[index - offset]
                offset += length
            raise Exception('index not found or not supported on superset')
    
        
    def _split_keys_over_sets(self, keys):
        split_sets = [ [] for x in self._private.particle_sets ]
        split_indices = [ [] for x in self._private.particle_sets ]
        
        
        if keys is None:
            offset = 0
            
            for setindex, x in enumerate(self._private.particle_sets):
                split_sets[setindex].extend(x.get_all_keys_in_store())
                split_indices[setindex].extend(range(offset, offset + len(x)))
                offset = offset + len(x)
                
        else:  
            if isinstance(keys, set):
                keys_array = numpy.array(list(keys))
            else:
                keys_array = numpy.array(keys)
                
            indices_array = numpy.arange(len(keys_array))
            for setindex, x in enumerate(self._private.particle_sets):
                mask = self._in1d(keys_array, x.get_all_keys_in_store(), True)
                split_sets[setindex] =  keys_array[mask]
                split_indices[setindex] = indices_array[mask]
        
        return split_sets, split_indices
                    
        
    def add_particles_to_store(self, keys, attributes = [], values = []):
        if not self._private.index_to_default_set is None:
            self._private.particle_sets[self._private.index_to_default_set].add_particles_to_store(keys, 
                attributes, values)
        else:
            raise exceptions.AmuseException("Cannot add particles to a superset")
        
    def remove_particles_from_store(self, keys):
        split_sets, split_indices = self._split_keys_over_sets(keys)
        for split_keys, set in zip(split_sets, self._private.particle_sets):
            set.remove_particles_from_store(split_keys)
    
    def get_values_in_store(self, keys, attributes):
        split_sets, split_indices = self._split_keys_over_sets(keys)
        
        indices_and_values = []
        
        for keys_for_set, indices_for_set, set in zip(split_sets, split_indices, self._private.particle_sets):
            if len(keys_for_set) > 0:
                values_for_set = set.get_values_in_store(keys_for_set, attributes)
                indices_and_values.append( (indices_for_set, values_for_set) )
        
        if keys is None:
            resultlength = len(self)
        else:
            resultlength = len(keys)
            
        values = [None] * len(attributes)
        units = [None] * len(attributes)
        for indices, values_for_set in indices_and_values:
            for valueindex, quantity in enumerate(values_for_set):
                resultvalue = values[valueindex]
                if resultvalue is None:
                    resultvalue = numpy.zeros(resultlength,dtype=quantity.number.dtype)
                    values[valueindex] = resultvalue
                    units[valueindex] = quantity.unit
                    
                resultunit = units[valueindex]
                
                numpy.put(resultvalue, indices, quantity.value_in(resultunit))
            
        return map(lambda u,v : u.new_quantity(v), units, values)
        
    def set_values_in_store(self, keys, attributes, values):
        split_sets, split_indices = self._split_keys_over_sets(keys)
        
        for keys_for_set, indices_for_set, set in zip(split_sets, split_indices, self._private.particle_sets):
            quantities = [None] * len(attributes)
            
            for valueindex, quantity in enumerate(values):
                if quantity.is_scalar():
                    numbers = [quantity.number]*len(indices_for_set)
                elif quantity.is_vector() and len(quantity) < len(keys):
                    numbers = numpy.take([quantity.number]*len(keys),indices_for_set)
                else:
                    numbers = numpy.take(quantity.number, indices_for_set)
                quantities[valueindex] = quantity.unit.new_quantity(numbers)
            
            set.set_values_in_store(keys_for_set, attributes, quantities)
    
    def get_attribute_names_defined_in_store(self):
        result = set(self._private.particle_sets[0].get_attribute_names_defined_in_store())
        for particle_set in self._private.particle_sets[1:]:
            result &= set(particle_set.get_attribute_names_defined_in_store())
        return list(result)
        
    def get_all_keys_in_store(self):
        result = []
        
        for set in self._private.particle_sets:
            result.extend(set.get_all_keys_in_store())
            
        return result
        
    def has_key_in_store(self, key):
        for set in self._private.particle_sets:
            if set.has_key_in_store(key):
                return True
        return False
        
        
    def _original_set(self):
        return self

    def _get_subset_for_key(self, key):
        for set in self._private.particle_sets:
            if set.has_key_in_store(key):
                return set
        return None

    def _is_superset(self):
        return True

    def _in1d(self, ar1, ar2, assume_unique=False):
        """
        copied from numpy.in1d (nump 1.4), to be compatible with numpy 1.3.0
        """
        if not assume_unique:
            ar1, rev_idx = numpy.unique(ar1, return_inverse=True)
            ar2 = numpy.unique(ar2)
    
        ar = numpy.concatenate( (ar1, ar2) )
        order = ar.argsort(kind='mergesort')
        sar = ar[order]
        equal_adj = (sar[1:] == sar[:-1])
        flag = numpy.concatenate( (equal_adj, [False] ) )
        index = order.argsort(kind='mergesort')[:len( ar1 )]
    
        if assume_unique:
            return flag[index]
        else:
            return flag[index][rev_idx]
            
    
    
    
class ParticlesSubset(AbstractParticleSet):
    """A subset of particles. Attribute values are not
    stored by the subset. The subset provides a limited view
    to the particles. 
    
    Particle subset objects are not supposed to be created
    directly. Instead use the ``to_set`` or ``select`` methods.

    """
    
    def __init__(self, particles, keys):
        AbstractParticleSet.__init__(self, particles)
        
        self._private.particles = particles
        self._private.keys = numpy.array(keys, dtype='uint64')
        self._private.set_of_keys = set(keys)
              
    
    def unconverted_set(self):
        return ParticlesSubset(self._private.particles.unconverted_set(), self._private.keys)
        
        
    def __getitem__(self, index):
        keys = self.get_all_keys_in_store()[index]
        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            key = keys
            if key == 0 or key >= (2**64 - 1):
                return None
            else:
                return Particle(self.get_all_keys_in_store()[index], self._original_set())
            
    def add_particles_to_store(self, keys, attributes = [], values = []):
        """
        Adds particles from to the subset, also
        adds the particles to the superset
        """
        self._private.keys = numpy.concatenate((self.keys,  numpy.array(keys,dtype='uint64')))
        self._private.set_of_keys += set(keys)
        self._private.particles.add_particles_to_store(keys, attributes, values)
        
    def remove_particles_from_store(self, keys):
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
    
    def get_values_in_store(self, keys, attributes):
        if keys is None:
            keys = self.get_all_keys_in_store()
            
        return self._private.particles.get_values_in_store(keys, attributes)
        
    def set_values_in_store(self, keys, attributes, values):
        if keys is None:
            keys = self.get_all_keys_in_store()
            
        self._private.particles.set_values_in_store(keys, attributes, values)
    
    def get_attribute_names_defined_in_store(self):
        return self._private.particles.get_attribute_names_defined_in_store()
        
    def get_all_keys_in_store(self):
        return self._private.keys
        
    def has_key_in_store(self, key):
        return key in self._private.set_of_keys        
        
    def _original_set(self):
        return self._private.particles
        
    def get_timestamp(self):
        return self._original_set().get_timestamp()
    
    def previous_state(self):
        return ParticlesSubset(self._private.particles.previous_state(), self._private.keys)
        
    def difference(self, other):
        new_set_of_keys = self._private.set_of_keys.difference(other.as_set()._private.set_of_keys)
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
        
        new_set_of_keys = self._private.set_of_keys.union(other.as_set()._private.set_of_keys)
        return ParticlesSubset(self._private.particles, list(new_set_of_keys))
    
    def as_set(self):
        return self



class ParticlesWithUnitsConverted(AbstractParticleSet):
    """
    A view on a particle sets. Used when to convert
    values between incompatible sets of units. For example to
    convert from si units to nbody units.
    
    The converter must have implement the ConverterInterface.
    
    >>> from amuse.units import nbody_system
    >>> particles_nbody = Particles(2)
    >>> particles_nbody.x = [10.0 , 20.0] | nbody_system.length
    >>> convert_nbody = nbody_system.nbody_to_si(10 | units.kg , 5 | units.m )
    >>> particles_si = ParticlesWithUnitsConverted(
    ...     particles_nbody,
    ...     convert_nbody.as_converter_from_si_to_nbody())
    ...
    >>> print particles_nbody.x
    [10.0, 20.0] length
    >>> print particles_si.x
    [50.0, 100.0] m
    >>> particles_si.x = [200.0, 400.0] | units.m
    >>> print particles_nbody.x
    [40.0, 80.0] length

    
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
        AbstractParticleSet.__init__(self, particles)
        
        self._private.particles = particles
        self._private.converter = converter
              
    def copy(self):
        copiedParticles =  self._private.particles.copy()
        return ParticlesWithUnitsConverted(copiedParticles, self._private.converter)
    
    def unconverted_set(self):
        return self._private.particles
        
    def add_particles_to_store(self, keys, attributes = [], values = []):
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_source_to_target(quantity)
            converted_values.append(converted_quantity)
        self._private.particles.add_particles_to_store(keys, attributes, converted_values)
        
    def remove_particles_from_store(self, keys):
        self._private.particles.remove_particles_from_store(keys)
        
    def get_values_in_store(self, keys, attributes):
        values = self._private.particles.get_values_in_store(keys, attributes)
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_target_to_source(quantity)
            converted_values.append(converted_quantity)
        return converted_values
        
    def set_values_in_store(self, keys, attributes, values):
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_source_to_target(quantity)
            converted_values.append(converted_quantity)
        self._private.particles.set_values_in_store(keys, attributes, converted_values)
    
    def get_attribute_names_defined_in_store(self):
        return self._private.particles.get_attribute_names_defined_in_store()
        
    def get_all_keys_in_store(self):
        return self._private.particles.get_all_keys_in_store()
        
    def has_key_in_store(self, key):
        return self._private.particles.has_key_in_store(key)
        
    def as_set(self):
        return ParticlesSubset(self, self.get_all_keys_in_store())
    
    def get_timestamp(self):
        timestamp = self._private.particles.get_timestamp()
        if not timestamp is None:
            timestamp = self._private.converter.from_target_to_source(timestamp)
        return timestamp
    
    def savepoint(self, timestamp=None):
        if not timestamp is None:
            timestamp = self._private.converter.from_target_to_source(timestamp)
        return ParticlesWithUnitsConverted(
            self._private.particles.savepoint(timestamp),
            self._private.converter
        )
    
    def previous_state(self):
        return ParticlesWithUnitsConverted(
            self._private.particles.previous_state(),
            self._private.converter
        )
    


    
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
        from_keys = self.from_particles.get_all_keys_in_store()
        to_keys = self.to_particles.get_all_keys_in_store()
        return numpy.intersect1d(from_keys,to_keys) #filter(lambda x : self.to_particles._has_key(x), from_keys)
        
    def copy_attributes(self, attributes):
        self._reindex()
        data = self.from_particles.get_values_in_store(self.keys, attributes)
        self.to_particles.set_values_in_store(self.keys, attributes, data)
    
    def copy(self):
        if not self.to_particles.can_extend_attributes():
            self.copy_overlapping_attributes()
        else:
            self.copy_all_attributes()
        
    def copy_all_attributes(self):
        names_to_copy = self.from_particles.get_attribute_names_defined_in_store()
        self._reindex()
        self.copy_attributes(list(names_to_copy))    

    def copy_overlapping_attributes(self):
        self._reindex()
        
        from_names = self.from_particles.get_attribute_names_defined_in_store()
        to_names = self.to_particles.get_attribute_names_defined_in_store()
        names_to_copy = set(from_names).intersection(set(to_names))
             
        self.copy_attributes(list(names_to_copy))  
        

    def copy_attribute(self, name, target_name = None):
        """ Copy the values of one attribute from the source set to the target set.
        The copied attribute can have a different name in the target set.
        
        :argument name: name of the attribute in the source set
        :argument target_name: name of the attribute in the target set, when None (default) the target_name
           will be set equal to the name
        
        >>> from amuse.datamodel import Particles
        >>> from amuse.units import units
        >>> particles1 = Particles(2)
        >>> particles2 = particles1.copy()
        >>> particles1.mass = 1 | units.m
        >>> particles2.mass = 3 | units.m
        >>> channel = particles1.new_channel_to(particles2)
        >>> channel.copy_attribute("mass", "mass_from_p2")
        >>> print particles2.mass_from_p2
        [1.0, 1.0] m
        >>> print particles2.mass - particles2.mass_from_p2
        [2.0, 2.0] m
        
        """
        if target_name is None:
            target_name = name
            
        self._reindex()
        data = self.from_particles.get_values_in_store(self.keys, [name,])
        self.to_particles.set_values_in_store(self.keys, [target_name,], data)
    
    
class Stars(Particles):

    def __init__(self, size = 0):
        Particles.__init__(self, size)

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
            if key == None:
                particles_set = Particles(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = Particles(1, keys = [key])
                
        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
       
        if isinstance(new_value_for_the_attribute, quantities.Quantity):
            self.particles_set._set_value_of_attribute(self.key, name_of_the_attribute, new_value_for_the_attribute)
        elif isinstance(new_value_for_the_attribute, Particle):
            self.particles_set._set_value_of_attribute(
                self.key, 
                name_of_the_attribute, 
                new_value_for_the_attribute.key | units.object_key
            )
        else:
            raise AttributeError("Can only assign quantities or other particles to an attribute.")
    
    def __getattr__(self, name_of_the_attribute):
        if hasattr(self.particles_set._private, 'particle_sets'):
            for current_set in self.particles_set._private.particle_sets:
                if current_set.has_key_in_store(self.key):
                    break
        else:
            current_set = self.particles_set
        try:
            return current_set._get_value_of_attribute(self.key, name_of_the_attribute)
        except Exception as ex:
            if ((name_of_the_attribute in current_set.get_attribute_names_defined_in_store()) or 
            (name_of_the_attribute in current_set._derived_attributes)):
                raise
            else:
                raise AttributeError("You tried to access attribute '{0}' but this attribute is not defined for this set.".format(name_of_the_attribute))
    
    def children(self):
        return self.particles_set.select(lambda x : x == self, ["parent"])
    
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
            raise exceptions.AmuseException("The parent and child particles should be in the same set")
        
        child.parent = self
                
    def __add__(self, particles):
        """
        Returns a particle subset, composed of the given
        particle(s) and this particle. Attribute values are
        not stored by the subset. The subset provides a view
        on the particles.
        
        :parameter particles: particle(s) to be added to self.
        
        >>> particles = Particles(2)
        >>> particle1 = particles[0]
        >>> particle1.x = 1.0 | units.m
        >>> particle2 = particles[1]
        >>> particle2.x = 2.0 | units.m
        >>> new_set = particle1 + particle2
        >>> new_set  # doctest:+ELLIPSIS
        <amuse.datamodel.particles.ParticlesSubset object at 0x...>
        >>> print len(new_set)
        2
        >>> print new_set.x
        [1.0, 2.0] m
        """
        return self.as_set().__add__(particles)
    
    def __sub__(self, particles):
        """
        Raises an exception: cannot subtract particle(s) 
        from a particle.
        """
        raise exceptions.AmuseException("Cannot subtract particle(s) from a particle.")
    
    def __str__(self):
        """
        Display string for a particle
        
        >>> p = Particle(10)
        >>> p.x = 10.2 | units.m
        >>> p.mass = 5 | units.kg
        >>> print p
        Particle(10, mass=5.0 kg, x=10.2 m)
        """
        output = 'Particle('
        output += str(self.key)
        for name, value in self.particles_set._values_of_particle(self.key):
            output += ', '
            output += name
            output += '='
            output += str(value)
        output += ')'
        return output
        
    
    def __dir__(self):
        result = []
        result.extend(dir(type(self)))
        result.extend(self.particles_set._attributes_for_dir())
        return result
        
    def __eq__(self, other):
        return isinstance(other, type(self)) and other.key == self.key
    
    def __hash__(self):
        return self.key.__hash__()
        
    def __ne__(self, other):
        return not (isinstance(other, type(self)) and other.key == self.key)

    def set_default(self, attribute, quantity):
        if not attribute in self.particles_set.get_attribute_names_defined_in_store():
            self.particles_set._set_value_of_attribute(self, attribute, quantity)
            
    def get_timeline_of_attribute(self, attribute):
        return self.particles_set.get_timeline_of_attribute(self.key, attribute)

    def get_timeline_of_attribute_as_vector(self, attribute):
        return self.particles_set.get_timeline_of_attribute_as_vector(self.key, attribute)

    def get_timeline_of_attributes(self, attributes):
        return self.particles_set.get_timeline_of_attributes(self.key, attributes)
        
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
    
    def as_particle_in_set(self, other):
        return other._get_particle(self.key)

def create_particle_set(**args):
    """
    Returns a particle set from the input vector quantities. input should be named
    keyword arguments.
    >>> m=units.kg([ 1.,1.])
    >>> x=units.m([0.,1.])
    >>> particles=create_particle_set(mass=m,x=x)
    >>> print len(particles)
    2
    """
    if len(args)==0:
        raise Exception("provide quantities")

    n=len(args.values()[0])

    for a in args:
        nn=len(args[a])
        if nn!=n:
            raise Exception("unequal length quantities")

    particles=Particles(n)

    for a in args:
        setattr(particles,a,args[a])

    return particles

