from amuse.support.core import CompositeDictionary
from amuse.support.core import compare_version_strings
from amuse.support import exceptions
from amuse.datamodel.base import *
from amuse.datamodel.memory_storage import *
from amuse.datamodel import trees
from amuse.units import constants
from amuse.units import units
from amuse.units import quantities
from amuse.units.quantities import Quantity
from amuse.units.quantities import new_quantity
from amuse.units.quantities import is_quantity
from amuse.units.quantities import as_vector_quantity
from amuse.units.quantities import zero
from amuse.units.quantities import AdaptingVectorQuantity

import random
import numpy
from numpy import ma


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
    quantity<[[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]] m>
    >>> particles.p[0]
    quantity<[1.0, 3.0, 5.0] m>
    >>> particles.position # "position" is a global vector attribute, coupled to x,y,z
    quantity<[[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]] m>

    """

    # this construct is needed to ensure that numpy see's grids
    # as objects and not as sequences
    # if we put a grid in a numpy object array we want the
    # grid in a field of that array and not the contents of
    # the grid (i.e. the grid points)
    # grids have the same trick
    if compare_version_strings(numpy.__version__, '1.7.0') < 0:
        __array_interface__ = {'shape':() }
    else:
        __array_interface__ = {'shape':(),'typestr':'|O4' }

    GLOBAL_DERIVED_ATTRIBUTES = {}


    def __init__(self, original = None):
        AbstractSet.__init__(self, original)



    #
    # Particle storage interface
    #

    def remove_particles_from_store(self, indices):
        pass

    def get_values_in_store(self, keys, attributes):
        pass

    def get_attribute_names_defined_in_store(self):
        return []

    def get_indices_of_keys(self, keys):
        pass


    #
    #
    #

    def _values_of_particle(self, index):
        attributes = self.get_attribute_names_defined_in_store()
        values = self.get_values_in_store(numpy.asarray([index]), attributes)

        for attribute, val in zip(attributes, values):
            yield attribute, val[0]

    #
    # public API
    #
    def __iter__(self):
        original_set = self._original_set()
        for key, index in zip(self.get_all_keys_in_store(), self.get_all_indices_in_store()):
            yield original_set._get_particle_unsave(key, index)


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


    def to_string(self, attributes_to_show = None, split_at = 20):
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
        >>> print particles.to_string()
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

        all_values = self.get_values_in_store(self.get_all_indices_in_store(), attributes)
        for index, quantity in enumerate(all_values):
            column = columns[index + 1]
            if hasattr(quantity, 'unit'):
                column.append(format_str11(str(quantity.unit)))
                quantity = quantity.number
            else:
                column.append(format_str11('none'))

            column.append('=' * 11)
            if len(quantity) > split_at * 2:
                if isinstance(quantity, LinkedArray):
                    values_to_show = list(map(format_str11,quantity[:split_at].to_print_list()))
                    values_to_show.append(format_str11('...'))
                    values_to_show.extend(map(format_str11,quantity[-split_at:].to_print_list()))
                elif hasattr(quantity, 'dtype'):
                    if numpy.issubdtype(quantity.dtype, float):
                        values_to_show = list(map(format_float,quantity[:split_at]))
                        values_to_show.append(format_str11('...'))
                        values_to_show.extend(map(format_float,quantity[-split_at:]))
                    else:
                        values_to_show = list(map(format_str11,quantity[:split_at]))
                        values_to_show.append(format_str11('...'))
                        values_to_show.extend(map(format_str11,quantity[-split_at:]))
                else:
                    values_to_show = list(map(format_str11,quantity[:split_at]))
                    values_to_show.append(format_str11('...'))
                    values_to_show.extend(map(format_str11,quantity[-split_at:]))
            else:
                if isinstance(quantity, LinkedArray):
                    values_to_show = map(format_str11,quantity.to_print_list())
                elif hasattr(quantity, 'dtype'):
                    if numpy.issubdtype(quantity.dtype, float):
                        try:
                            values_to_show = map(format_float,quantity)
                        except ValueError:
                            values_to_show = map(format_str11,quantity)
                    else:
                        values_to_show = map(format_str11,quantity)
                else:
                    values_to_show = map(format_str11, quantity)

            column.extend(values_to_show)
            column.append('=' * 11)

        column = columns[0]
        column.append(format_str20("-"))
        column.append('=' * 20)
        particle_keys = self.get_all_keys_in_store()
        if len(particle_keys) > split_at * 2:
            values_to_show = list(map(format_str20, particle_keys[:split_at]))
            values_to_show.append(format_str20('...'))
            values_to_show.extend(map(format_str20, particle_keys[-split_at:]))
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

    def _get_particle_unsave(self, key, index = -1):
        return Particle(
            key,
            self._original_set(),
            set_index = index,
            set_version = self._get_version()
        )

    def can_extend_attributes(self):
        return self._original_set().can_extend_attributes()

    def add_attribute_domain(self, namespace):
        self._derived_attributes[namespace] = DomainAttribute(namespace)

    def _is_superset(self):
        return False

    def as_binary_tree(self, name_of_first_child = 'child1', name_of_second_child = 'child2'):
        return trees.ChildTreeOnParticleSet(self, (name_of_first_child, name_of_second_child))

    def new_binary_tree_wrapper(self, name_of_first_child = 'child1', name_of_second_child = 'child2'):
        return trees.ChildTreeOnParticleSet(self, (name_of_first_child, name_of_second_child))

    def copy(self, memento = None, keep_structure = False, filter_attributes = lambda particle_set, x : True):
        attributes = self.get_attribute_names_defined_in_store()
        attributes = [x for x in attributes if filter_attributes(self, x)]
        keys = self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        values = self.get_values_in_store(indices, attributes)
        result = self._factory_for_new_collection()()
        if memento is None:
            memento = {}
        memento[id(self._original_set())] = result

        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy(memento, keep_structure, filter_attributes))
            else:
                converted.append(x)
        result.add_particles_to_store(keys, attributes, converted)

        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))

        return result


    def copy_to_new_particles(self, keys = None, keys_generator = None, memento = None, keep_structure = False, filter_attributes = lambda particle_set, x : True):
        if keys_generator is None:
            keys_generator = UniqueKeyGenerator

        my_keys = self.get_all_keys_in_store()

        if not keys is None:
            if len(keys) != len(my_keys):
                    raise Exception('not enough new keys given for the copy')
            else:
                particle_keys = keys
        else:
            particle_keys = keys_generator.next_set_of_keys(len(my_keys))

        attributes = self.get_attribute_names_defined_in_store()
        attributes = [x for x in attributes if filter_attributes(self, x)]
        indices = self.get_all_indices_in_store()
        values = self.get_values_in_store(indices, attributes)
        result = self._factory_for_new_collection()()
        if memento is None:
            memento = {}

        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy(memento, keep_structure, filter_attributes))
            else:
                converted.append(x)

        memento[id(self._original_set())] = result

        result.add_particles_to_store(particle_keys, attributes, converted)

        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))

        return result

    def _factory_for_new_collection(self):
        return Particles

    def empty_copy(self):
        """
        Creates a new in memory set and copies the particles to it.
        The attributes and values are not copied.The history
        of the set is not copied over.

        >>> from amuse.datamodel import Particles
        >>> from amuse.units import units
        >>> original = Particles(2)
        >>> original.mass = 0 | units.m
        >>> print hasattr(original, "mass")
        True
        >>> print len(original)
        2
        >>> copy = original.empty_copy()
        >>> print hasattr(copy, "mass")
        False
        >>> print len(copy)
        2

        """
        keys = self.get_all_keys_in_store()
        result = Particles()
        result.add_particles_to_store(keys, [],[])
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

    def copy_values_of_attributes_to(self, attribute_names, particles):
        channel = self.new_channel_to(particles)
        channel.copy_attributes(attribute_names)

    def new_channel_to(self, other, attributes=None, target_names=None):
        return ParticleInformationChannel(self, other, attributes, target_names)

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
        if not original_particles_set is particles._original_set():
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
                raise exceptions.AmuseException("Unable to subtract a particle, because it is not part of this set.")
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
        indices = particles.get_all_indices_in_store()
        keys =  particles.get_all_keys_in_store()
        values = particles.get_values_in_store(indices, attributes)
        values = map(self._convert_from_entities_or_quantities, values)
        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy_with_link_transfer(particles, self))
            else:
                converted.append(x)
        try:
            self.add_particles_to_store(keys, attributes, converted)
        except exceptions.MissingAttributesAmuseException as caught_exception:
            for attribute_name in caught_exception.missing_attributes:
                if attribute_name in particles._derived_attributes:
                    attributes.append(attribute_name)
                    converted.append(getattr(particles, attribute_name))
                else:
                    raise
            self.add_particles_to_store(keys, attributes, converted)
        return ParticlesSubset(self._original_set(),keys)


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
        indices = self.get_indices_of_keys(particles.get_all_keys_in_store())
        self.remove_particles_from_store(indices)


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
            values = self.get_values_in_store(self.get_indices_of_keys(added_keys), attributes)
            converted = []
            for x in values:
                if isinstance(x, LinkedArray):
                    converted.append(x.copy_with_link_transfer(self._original_set(), other_particles))
                else:
                    converted.append(x)
            other_particles.add_particles_to_store(added_keys, attributes, converted)

        removed_keys = list(removed_keys)
        if removed_keys:
            other_particles.remove_particles_from_store(other_particles.get_indices_of_keys(removed_keys))

    def compressed(self):
        return self

    def get_valid_particles_mask(self):
        return numpy.ones(len(self), dtype = numpy.bool)

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

    def _masked_subset(self, keys):
        return ParticlesMaskedSubset(self._original_set(), keys)

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
        return self.__getitem__(getattr(self, attribute).argsort(kind=kind))


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
        indices = self.get_all_indices_in_store()
        keys = self.get_all_keys_in_store()
        values = self.get_values_in_store(indices, attributes)
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


    def __contains__(self, particle):
        """
        Check if a particle is part of a set

        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> p1 = particles[1]
        >>> print p1 in particles
        True
        >>> p4 = Particle()
        >>> print p4 in particles
        False

        """
        keys = set(self.get_all_keys_in_store())
        return particle.key in keys


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


    def _as_masked_subset_in(self, other):
        keys = numpy.ma.array(self.get_all_keys_in_store(), dtype='uint64')
        keys.mask = ~self.get_valid_particles_mask()
        return other._masked_subset(keys)

    def random_sample(self, number_of_particles):
        return self.__getitem__(random.sample(xrange(len(self)), number_of_particles))


class Particles(AbstractParticleSet):
    """
    A set of particles. Attributes and values are stored in
    a private storage model. This storage model can store
    the values in the python memory space, in the memory space
    of the code or in a HDF5 file. By default the storage
    model is in memory.



    """
    def __init__(self, size = 0, storage = None, keys = None, keys_generator = None, particles = None, is_working_copy = True, **attributes):
        AbstractParticleSet.__init__(self)

        self._private.version = 0
        self._private.is_working_copy = is_working_copy

        if storage is None:
            self._private.attribute_storage = get_in_memory_attribute_storage_factory()()
        else:
            self._private.attribute_storage = storage


        if keys_generator is None:
            keys_generator = UniqueKeyGenerator

        if not particles is None:
            if isinstance(particles,AbstractParticleSet):
                self.add_particles(particles)
            else:
                for x in iter(particles):
                    self.add_particle(x)
        elif size > 0:
            if not keys is None:
                if len(keys) != size:
                    raise Exception('keys and size was specified in the creation of a particle set, but the length of the keys is not equal to the size')
                else:
                    particle_keys = keys
            else:
                particle_keys = keys_generator.next_set_of_keys(size)

            self.add_particles_to_store(particle_keys)
        elif not keys is None:
            self.add_particles_to_store(keys)
        elif len(attributes) > 0:
            number_of_attributes = 0
            for attributevalue in attributes.values():
                if is_quantity(attributevalue):
                    if attributevalue.is_scalar():
                        number_of_attributes = max(number_of_attributes, 1)
                    else:
                        number_of_attributes = max(number_of_attributes, len(attributevalue))
                else:
                    try:
                        if isinstance(attributevalue, basestring):
                            number_of_attributes = max(number_of_attributes,1)
                        else:
                            number_of_attributes = max(number_of_attributes, len(attributevalue))
                    except: #fails for numbers
                        number_of_attributes = max(number_of_attributes,1)

            particle_keys = keys_generator.next_set_of_keys(number_of_attributes)
            self.add_particles_to_store(particle_keys)


        self._private.version = 0

        if len(attributes) > 0:
            attributenames = []
            attributevalues = []
            for attributename, attributevalue in attributes.iteritems():
                attributenames.append(attributename)
                attributevalues.append(attributevalue)

            self.set_values_in_store(self.get_all_indices_in_store(), attributenames, attributevalues)

        self._private.previous = None
        self.collection_attributes.timestamp = None


    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]
        index = self.get_all_indices_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Particle(keys, self, index, self._get_version())

    def _get_version(self):
        return self._private.version

    def __iter__(self):
        keys =  self.get_all_keys_in_store()
        indices = self.get_all_indices_in_store()
        version = self._get_version()

        for i in range(len(keys)):
            yield Particle(keys[i], self,  indices[i], version)

    def savepoint(self, timestamp=None, format = 'memory', **attributes):
        if format == 'memory':
            instance = type(self)(is_working_copy = False)
            instance._private.attribute_storage = self._private.attribute_storage.copy()
            instance.collection_attributes.timestamp = timestamp

            for name, value in attributes.iteritems():
                setattr(instance.collection_attributes, name, value)
        else:
            raise Exception("{0} not supported, only 'memory' savepoint supported".format(format))
        instance._private.previous = self._private.previous
        instance._private.version = 0
        self._private.previous = instance
        return instance

    def new_working_copy(self):
        if self._private.is_working_copy:
            previous = self._private.previous
            if previous is None:
                raise Exception("you have not savepoint for this set, you cannot create a working copy please use copy instead".format(format))
        else:
            previous = self
        result = previous.copy()
        result._private.previous = previous
        return result

    def get_timestamp(self):
        return self.collection_attributes.timestamp

    def iter_history(self):
        if self._private.is_working_copy:
            current = self._private.previous
        else:
            current = self

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
                index = x.get_indices_of_keys([particle_key])[0]
                timeline.append((x.collection_attributes.timestamp, x._get_value_of_attribute(x[index], index, attribute)))
        return timeline

    def get_timeline_of_attribute_as_vector(self, particle_key, attribute):
        timeline = AdaptingVectorQuantity()
        chrono_values = AdaptingVectorQuantity()
        for x in self.history:
            if x.has_key_in_store(particle_key):
                index = x.get_indices_of_keys([particle_key])[0]
                timeline.append(x.collection_attributes.timestamp)
                chrono_values.append(x._get_value_of_attribute(x[index], index, attribute))
        return timeline, chrono_values

    def get_timeline_of_attributes(self, particle_key, attributes):
        result = map(lambda x: [], range(len(attributes)+1))
        units = map(lambda x: None, range(len(attributes)+1))

        for x in self.history:
            if x.has_key_in_store(particle_key):
                index = x.get_indices_of_keys([particle_key])[0]
                if  units[0] is None:
                    units[0] = x.collection_attributes.timestamp.unit
                for i, attribute in enumerate(attributes):
                    quantity = x._get_value_of_attribute(x[index], index, attribute)
                    if  units[i+1] is None:
                        units[i+1] = quantity.unit
                    result[i+1].append(quantity.value_in(units[i+1]))

        return list(map(lambda value,unit : unit.new_quantity(value), result, units))


    def add_particles_to_store(self, keys, attributes = [], values = []):
        self._private.attribute_storage.add_particles_to_store(keys, attributes, values)
        self._private.version += 1


    def remove_particles_from_store(self, indices):
        self._private.attribute_storage.remove_particles_from_store(indices)
        self._private.version += 1

    def get_values_in_store(self, indices, attributes):
        missing_attributes = set(attributes) - set(self.get_attribute_names_defined_in_store()) - set(["index_in_code"])

        if len(missing_attributes) == 0:
            return self._private.attribute_storage.get_values_in_store(indices, attributes)

        defined_attributes = list(set(attributes) - missing_attributes)
        defined_values = dict(zip(
            defined_attributes,
            self._private.attribute_storage.get_values_in_store(indices, defined_attributes)
        ))
        subset = self[indices]
        return [defined_values[attribute] if attribute in defined_values else subset._get_derived_attribute_value(attribute) for attribute in attributes]

    def get_indices_of_keys(self, keys):
        return self._private.attribute_storage.get_indices_of(keys)

    def set_values_in_store(self, indices, attributes, values):
        self._private.attribute_storage.set_values_in_store(indices, attributes, values)

    def get_attribute_names_defined_in_store(self):
        return self._private.attribute_storage.get_defined_attribute_names()

    def get_all_keys_in_store(self):
        return self._private.attribute_storage.get_all_keys_in_store()

    def get_all_indices_in_store(self):
        return self._private.attribute_storage.get_all_indices_in_store()

    def has_key_in_store(self, key):
        return self._private.attribute_storage.has_key_in_store(key)

    def get_value_in_store(self, index, attribute):
        return self._private.attribute_storage.get_value_in_store(index, attribute)

    def can_extend_attributes(self):
        return self._private.attribute_storage.can_extend_attributes()


    def _remove_indices_in_attribute_storage(self, indices):
        self._private.attribute_storage._remove_indices(indices)
        self._private.version += 1

    def _add_indices_in_attribute_storage(self, indices):
        self._private.attribute_storage._add_indices(indices)
        self._private.version += 1

class BoundSupersetParticlesFunctionAttribute(object):
    def  __init__(self, name, superset):
        self.name = name
        self.superset = superset
        self.subsetfunctions = []

    def add_subsetfunction(self, callable):
        self.subsetfunctions.append(callable)

    def __call__(self, *list_arguments, **keyword_arguments):
        subset_results = []
        for x in self.subsetfunctions:
            subset_results.append(x(*list_arguments, **keyword_arguments))

        if subset_results[0] is None:
            return None
        if isinstance(subset_results[0], AbstractParticleSet):
            return ParticlesSuperset(subset_results)
        if hasattr(subset_results[0], 'unit'):
            result = AdaptingVectorQuantity()
            for one_result in subset_results:
                result.extend(one_result)
            return result
        return [item for one_result in subset_results for item in one_result]

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

            elif hasattr(subset_result, 'dtype'):
                if len(subset_result) == 0:
                    continue
                if result is None:
                    shape = [len(superset),] + list(subset_result.shape[1:])
                    result = numpy.zeros(shape, dtype=subset_result.dtype)
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

    def get_value_for_entity(self, superset, particle, index):
        raise exceptions.AmuseException("Internal AMUSE error, a single entity (Particle) should always be bound to the subset and not the superset")

    def set_value_for_entity(self, superset, key, value):
        raise exceptions.AmuseException("Internal AMUSE error, a single entity (Particle) should always be bound to the subset and not the superset")

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

    def __init__(self, particle_sets, index_to_default_set=None, names = None):
        AbstractParticleSet.__init__(self)

        if not names is None:
            self._private.mapping_from_name_to_set = {}
            for name, particle_set in zip(names, particle_sets):
                self._private.mapping_from_name_to_set[name] = particle_set

        self._private.particle_sets = list(particle_sets)
        self._private.index_to_default_set = index_to_default_set


        names_of_derived_attributes_in_all_subsets = None
        for subset in particle_sets:
            derived_attribute_names = set(subset._derived_attributes.keys())
            if names_of_derived_attributes_in_all_subsets is None:
                names_of_derived_attributes_in_all_subsets = set(derived_attribute_names)
            else:
                names_of_derived_attributes_in_all_subsets &= derived_attribute_names

        names_of_derived_attributes_in_all_subsets -= set(self.GLOBAL_DERIVED_ATTRIBUTES.keys())
        for name in names_of_derived_attributes_in_all_subsets:
            self._derived_attributes[name] = DerivedSupersetAttribute(name)

        self._private.version = -1
        self._ensure_updated_set_properties()

        if self.has_duplicates():
            raise exceptions.AmuseException("Unable to add a particle, because it was already part of this set.")


    def _ensure_updated_set_properties(self):
        if self._private.version == self._get_subsets_version():
            return

        self._private.version = self._get_subsets_version()
        self._private.length = numpy.sum([len(x) for x in  self._private.particle_sets])
        self._private.indices = numpy.arange(self._private.length)
        self._private.keys = self._get_concatenated_keys_in_store()
        self._private.key_to_index = {}

        d = self._private.key_to_index
        index = 0
        for x in self._private.keys:
            d[x] = index
            index += 1

    def can_extend_attributes(self):
        for x in self._private.particle_sets:
            if not x.can_extend_attributes():
                return False
        return True

    def __len__(self):
        self._ensure_updated_set_properties()

        return self._private.length

    def __iter__(self):
        for set in self._private.particle_sets:
            for particle in set:
                yield particle

    def _get_subsets_version(self):
        versions = [[x._get_version()] for x in self._private.particle_sets]
        return numpy.sum(versions)

    def _get_version(self):
        self._ensure_updated_set_properties()

        return self._private.version


    def __getitem__(self, index):
        self._ensure_updated_set_properties()

        offset = 0

        if isinstance(index, basestring):
            return self.get_subset(index)
        else:
            keys = self.get_all_keys_in_store()[index]
            if hasattr(keys, '__iter__'):
                return self._subset(keys)
            else:
                index = self.get_all_indices_in_store()[index]
                for set in self._private.particle_sets:
                    length = len(set)
                    if index < (offset+length):
                        return set[index - offset]
                    offset += length
                raise Exception('index not found on superset')

    def _get_particle(self, key):
        if self.has_key_in_store(key):
            return self._get_subset_for_key(key)._get_particle(key)
        else:
            return None

    def _get_particle_unsave(self, key, index = -1):
        if index >= 0:
            offset, subset = self._get_offset_and_subset_for_index(index)
            index -= offset
            return subset._get_particle_unsave(key, subset.get_all_indices_in_store()[index])
        else:
            return self._get_subset_for_key(key)._get_particle_unsave(key)

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


    def _split_indices_over_sets(self, indices):
        self._ensure_updated_set_properties()

        split_sets = [ [] for x in self._private.particle_sets ]
        split_indices = [ [] for x in self._private.particle_sets ]

        offset = 0
        if isinstance(indices, set):
            indices = numpy.array(list(indices))

        if indices is None:
            offset = 0

            for setindex, x in enumerate(self._private.particle_sets):
                split_sets[setindex].extend(x.get_all_indices_in_store())
                split_indices[setindex].extend(range(offset, offset + len(x)))
                offset = offset + len(x)
        elif len(indices) == 0:
            for setindex, x in enumerate(self._private.particle_sets):
                split_sets[setindex] = []
                split_indices[setindex] = []
        else:
            result_indices_array = numpy.arange(len(indices))
            for setindex, x in enumerate(self._private.particle_sets):
                mask = numpy.logical_and( (indices >= offset) , (indices < (offset + len(x))) )
                indices_in_store = numpy.asarray(x.get_all_indices_in_store())
                split_sets[setindex] = indices_in_store[(indices-offset)[mask]]
                split_indices[setindex] = result_indices_array[mask]
                offset = offset + len(x)
        return split_sets, split_indices


    def add_particles_to_store(self, keys, attributes = [], values = []):
        if not self._private.index_to_default_set is None:
            self._private.particle_sets[self._private.index_to_default_set].add_particles_to_store(keys,
                attributes, values)
        else:
            raise exceptions.AmuseException("Cannot add particles to a superset")

    def remove_particles_from_store(self, indices):
        split_indices_in_subset, split_indices_in_input = self._split_indices_over_sets(indices)

        for indices_in_subset, set in zip(split_indices_in_subset, self._private.particle_sets):
            if len(indices_in_subset) == 0:
                continue

            set.remove_particles_from_store(indices_in_subset)


    def get_values_in_store(self, indices, attributes):
        split_indices_in_subset, split_indices_in_input = self._split_indices_over_sets(indices)

        indices_and_values = []

        for indices_in_subset, indices_in_input, set in zip(split_indices_in_subset, split_indices_in_input, self._private.particle_sets):
            if len(indices_in_subset) > 0:
                values_for_set = set.get_values_in_store(indices_in_subset, attributes)
                indices_and_values.append( (indices_in_input, values_for_set) )

        if indices is None:
            resultlength = len(self)
        else:
            resultlength = len(indices)

        values = [[]] * len(attributes)
        units = [None] * len(attributes)
        converts = [lambda x : x] * len(attributes)
        for indices, values_for_set in indices_and_values:
            for valueindex, quantity in enumerate(values_for_set):
                resultvalue = values[valueindex]
                if len(resultvalue) == 0:
                    if not is_quantity(quantity):
                        dtype = quantity.dtype
                        converts[valueindex] = lambda x : x
                        units[valueindex] = None
                    else:
                        dtype = quantity.number.dtype
                        converts[valueindex] = quantity.unit.new_quantity
                        units[valueindex] = quantity.unit
                    shape = list(quantity.shape)
                    shape[0] = resultlength
                    resultvalue = numpy.zeros(shape,dtype=dtype)
                    values[valueindex] = resultvalue

                resultunit = units[valueindex]
                if not resultunit is None:
                    resultvalue[indices] = quantity.value_in(resultunit)
                else:
                    resultvalue[indices] = quantity

        return map(lambda u,v : u(v), converts, values)

    def set_values_in_store(self, indices, attributes, values):
        split_indices_in_subset, split_indices_in_input = self._split_indices_over_sets(indices)

        for indices_in_subset, indices_in_input, set in zip(split_indices_in_subset, split_indices_in_input, self._private.particle_sets):
            quantities = [None] * len(attributes)

            for valueindex, quantity in enumerate(values):
                if is_quantity(quantity):
                    if quantity.is_scalar():
                        numbers = [quantity.number]*len(indices_in_input)
                    elif quantity.is_vector() and len(quantity) < len(indices):
                        numbers = numpy.take([quantity.number]*len(indices),indices_in_input)
                    else:
                        numbers = numpy.take(quantity.number, indices_in_input)
                    quantities[valueindex] = quantity.unit.new_quantity(numbers)
                else:
                    if not hasattr(quantity, 'ndim'):
                        numbers = numpy.asarray([quantity]*len(indices_in_input))
                    elif len(quantity) < len(indices):
                        numbers = numpy.take([quantity]*len(indices),indices_in_input)
                    else:
                        numbers = numpy.take(quantity, indices_in_input)
                    quantities[valueindex] = numbers

            set.set_values_in_store(indices_in_subset, attributes, quantities)

    def get_attribute_names_defined_in_store(self):
        self._ensure_updated_set_properties()
        result = set(self._private.particle_sets[0].get_attribute_names_defined_in_store())
        for particle_set in self._private.particle_sets[1:]:
            if len(particle_set) > 0:
                result &= set(particle_set.get_attribute_names_defined_in_store())
        return list(result)

    def get_all_keys_in_store(self):
        self._ensure_updated_set_properties()

        return self._private.keys

    def get_all_indices_in_store(self):
        self._ensure_updated_set_properties()

        return self._private.indices

    def get_indices_of_keys(self, keys):
        self._ensure_updated_set_properties()

        return numpy.array([self._private.key_to_index[x] for x in keys])


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

    def _get_offset_and_subset_for_index(self, index):
        offset = 0
        for set in self._private.particle_sets:
            setlen = len(set)
            if index >= offset and index < (offset+setlen):
                return offset, set
            offset += setlen
        return None, None

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

    def _get_concatenated_keys_in_store(self):
        result = []
        dtype = None
        for set in self._private.particle_sets:
            subsetkeys = set.get_all_keys_in_store()
            if dtype is None and len(set) > 0:
                dtype = subsetkeys.dtype
            result.extend(subsetkeys)

        return numpy.array(result, dtype = dtype)

    def get_subsets(self):
        return list(self._private.particle_sets)

    def get_subset(self, name):
        return self._private.mapping_from_name_to_set[name]

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

        self._private.version = -1
        self._private.indices = None

    def __getitem__(self, index):
        keys = self.get_all_keys_in_store()[index]
        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            key = keys
            if key > 0 and key < 18446744073709551615L: #2**64 - 1
                return self._private.particles._get_particle_unsave(key, self.get_all_indices_in_store()[index])
            else:
                return None

    def _get_version(self):
        return self._private.particles._get_version()

    def compressed(self):
        keys = self._private.keys
        return self._subset(keys[numpy.logical_and(keys > 0 ,  keys < 18446744073709551615L)])

    def get_valid_particles_mask(self):
        keys = self._private.keys
        return numpy.logical_and(keys > 0 ,  keys < 18446744073709551615L)

    def unconverted_set(self):
        return ParticlesSubset(self._private.particles.unconverted_set(), self._private.keys)

    def add_particles_to_store(self, keys, attributes = [], values = []):
        """
        Adds particles from to the subset, also
        adds the particles to the superset
        """
        self._private.keys = numpy.concatenate((self._private.keys,  numpy.array(keys,dtype='uint64')))
        self._private.set_of_keys |= set(keys)
        self._private.particles.add_particles_to_store(keys, attributes, values)

    def remove_particles_from_store(self, indices):
        """
        Removes particles from the subset, and removes particles
        from the super set
        """

        indices_to_remove = set(indices)

        index = 0
        index_in_local_list = []
        for x in self._private.indices:
            if x in indices_to_remove:
                index_in_local_list.append(index)
            index += 1
        set_to_remove = set(self._private.keys[index_in_local_list])
        self._private.keys = numpy.delete(self._private.keys,index_in_local_list)
        self._private.set_of_keys -= set_to_remove
        self._private.particles.remove_particles_from_store(indices)

        self._private.version = -1
        self._private.indices = None


    def get_values_in_store(self, indices, attributes):
        if indices is None:
            indices = self.get_all_indices_in_store()

        return self._private.particles.get_values_in_store(indices, attributes)

    def set_values_in_store(self, indices, attributes, values):
        if indices is None:
            indices = self.get_all_indices_in_store()

        self._private.particles.set_values_in_store(indices, attributes, values)

    def get_attribute_names_defined_in_store(self):
        return self._private.particles.get_attribute_names_defined_in_store()

    def get_all_keys_in_store(self):
        return self._private.keys

    def get_all_indices_in_store(self):
        if not self._private.version == self._private.particles._get_version():
            self._private.indices = self._private.particles.get_indices_of_keys(self._private.keys)
            self._private.version = self._private.particles._get_version()

        return self._private.indices

    def has_key_in_store(self, key):
        return key in self._private.set_of_keys

    def _original_set(self):
        return self._private.particles

    def get_timestamp(self):
        return self._original_set().get_timestamp()

    def get_indices_of_keys(self, keys):
        if not self._private.version == self._private.particles._get_version():
            self._private.indices = self._private.particles.get_indices_of_keys(self._private.keys)
            self._private.version = self._private.particles._get_version()

        return self._private.particles.get_indices_of_keys(keys)

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

    def copy(self, memento = None, keep_structure = False, filter_attributes = lambda particle_set, x : True):
        if keep_structure:
            result = ParticlesSubset(None, [])
            if memento is None:
                memento = dict()
            memento[id(self)] = result
            if id(self._private.particles) in memento:
                result._private.particles = memento[id(self._private.particles)]
            else:
                result._private.particles = self._private.particles.copy(memento, keep_structure, filter_attributes)
            result._private.keys = numpy.array(self._private.keys, dtype='uint64')
            result._private.set_of_keys = set(self._private.keys)

            result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
            object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))

            return result
        else:
            return super(ParticlesSubset, self).copy(memento, keep_structure, filter_attributes)



class ParticlesMaskedSubset(ParticlesSubset):
    """A subset of particles. Attribute values are not
    stored by the subset. The subset provides a limited view
    to the particles. Unlike the normal subset the masked subset
    can store None values.
    """

    def __init__(self, particles, keys):
        AbstractParticleSet.__init__(self, particles)

        self._private.particles = particles

        self._private.keys = numpy.ma.array(keys, dtype='uint64').copy()
        validkeys = self._private.keys.compressed()
        if len(validkeys) > 0:
            self._private.set_of_keys = set(validkeys)
        else:
            self._private.set_of_keys = set([])

        self._private.version = -1
        self._private.indices = None

    def compressed(self):
        return self._subset(self._private.keys.compressed())

    def get_valid_particles_mask(self):
        return ~self._private.keys.mask

    def __iter__(self):
        original_set = self._original_set()
        for key in self._private.keys :
            if key is ma.masked:
                yield None
            else:
                yield original_set._get_particle_unsave(key)

    def __getitem__(self, index):
        keys = self.get_all_keys_in_store()[index]

        if not keys is ma.masked and hasattr(keys, '__iter__'):
            if numpy.all(~ keys.mask):
                return self._subset(keys.data)
            else:
                return ParticlesMaskedSubset(self._original_set(),keys)
        else:
            key = keys
            if not key is ma.masked:
                return self._original_set()._get_particle_unsave(key, self.get_all_indices_in_store()[index])
            else:
                return None

    def _get_version(self):
        return self._private.particles._get_version()

    def unconverted_set(self):
        return ParticlesMaskedSubset(self._private.particles.unconverted_set(), self._private.keys)

    def add_particles_to_store(self, keys, attributes = [], values = []):
        raise Exception("Cannot add particles to a masked subset")

    def remove_particles_from_store(self, keys):
        raise Exception("Cannot remove particles from a masked subset")


    def get_values_in_store(self, indices, attributes):
        if indices is None:
            indices = self.get_all_indices_in_store()
        return self._private.particles.get_values_in_store(indices, attributes)


        #subkeys = keys[~mask]
        #subvalues = self._private.particles.get_values_in_store(subkeys, attributes)

        #raise Exception("not implemented more of this")


    def set_values_in_store(self, indices, attributes, values):
        if indices is None:
            indices = self.get_all_indices_in_store()
            mask = self._private.keys.mask
        else:
            mask = indices.mask

        if numpy.all(~mask):
            return self._private.particles.set_values_in_store(indices, attributes, values)

        raise Exception("not implemented more of this")

    def get_attribute_names_defined_in_store(self):
        return self._private.particles.get_attribute_names_defined_in_store()

    def get_all_keys_in_store(self):
        return self._private.keys

    def get_all_indices_in_store(self):
        if not self._private.version == self._private.particles._get_version():
            self._private.indices = -1 * numpy.ones(len(self._private.keys), dtype = numpy.int64)
            self._private.indices = numpy.ma.array(self._private.indices, dtype = numpy.int64)
            self._private.indices.mask = self._private.keys.mask
            mask = self._private.keys.mask
            mask = ~mask
            self._private.indices[mask] = self._private.particles.get_indices_of_keys(self._private.keys[mask])
            self._private.version = self._private.particles._get_version()

        return self._private.indices

    def has_key_in_store(self, key):
        return key in self._private.set_of_keys

    def get_indices_of_keys(self, keys):
        if len(keys) == 1:
            return self._private.particles.get_indices_of_keys(keys)
        keys = numpy.ma.array(keys, dtype='uint64')
        result = -1 * numpy.ones(len(keys), dtype = numpy.int32)
        result[~keys.mask] = self._private.particles.get_indices_of_keys(keys[~keys.mask])
        return result

    def previous_state(self):
        return ParticlesMaskedSubset(self._private.particles.previous_state(), self._private.keys)


    def copy(self, memento = None, keep_structure = False, filter_attributes = lambda particle_set, x : True):
        attributes = self.get_attribute_names_defined_in_store()
        attributes = [x for x in attributes if filter_attributes(self, x)]

        keys = self.get_all_keys_in_store()
        keys = keys[~keys.mask]
        values = self.get_values_in_store(keys, attributes)
        result = Particles()
        converted = []
        if memento is None:
            memento = {}
        memento[id(self._original_set())] = result
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy(memento, keep_structure, filter_attributes))
            else:
                converted.append(x)
        result.add_particles_to_store(keys, attributes, converted)
        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))
        return result



    def copy_to_new_particles(self, keys = None, keys_generator = None, memento = None, keep_structure = False, filter_attributes = lambda particle_set, x : True):
        if keys_generator is None:
            keys_generator = UniqueKeyGenerator

        my_keys = self.get_all_keys_in_store()
        my_keys = my_keys[~my_keys.mask]

        if not keys is None:
            if len(keys) != len(my_keys):
                raise Exception('not enough new keys given for the copy')
            else:
                particle_keys = keys
        else:
            particle_keys = keys_generator.next_set_of_keys(len(my_keys))

        attributes = self.get_attribute_names_defined_in_store()
        attributes = [x for x in attributes if filter_attributes(self, x)]
        values = self.get_values_in_store(my_keys, attributes)
        result = Particles()
        if memento is None:
            memento = {}

        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy(memento, keep_structure, filter_attributes))
            else:
                converted.append(x)

        memento[id(self._original_set())] = result

        result.add_particles_to_store(particle_keys, attributes, converted)

        result._private.collection_attributes = self._private.collection_attributes._copy_for_collection(result)
        object.__setattr__(result, "_derived_attributes", CompositeDictionary(self._derived_attributes))

        return result

    def __str__(self):
        """
        Display string of a particle subset.

        >>> particles = Particles(keys=[1,2])
        >>> particles[0].child = particles[1]
        >>> print particles.child[0]
        Particle(2, child=None)
        >>> print particles.child[1]
        None
        """
        keys = numpy.where(self.get_valid_particles_mask(), self._private.keys, [None])
        return str(list(keys))

    def as_set(self):
        return self


class ParticlesOverlay(AbstractParticleSet):
    """An overlay of of a particles set. The overlay
    stores extra attributes for the particles in the
    overlayed set.

    >>> p1 = Particles(3)
    >>> p1.mass = [10.0, 20.0, 30.0] | units.kg
    >>> p2 = ParticlesOverlay(p1)
    >>> p2.radius = [4.0, 5.0, 6.0] | units.m
    >>> print len(p2)
    3
    >>> print p2.mass
    [10.0, 20.0, 30.0] kg
    """

    def __init__(self, particles, overlay_set = None):
        AbstractParticleSet.__init__(self)

        self._private.base_set = particles
        if overlay_set is None:
            overlay_set = Particles(keys = self._private.base_set.key)
        self._private.overlay_set = overlay_set
        self._private.base_version = self._private.base_set._get_version()

    def _ensure_updated_set_properties(self):
        if self._private.base_version == self._private.base_set._get_version():
            return

        self._private.base_version = self._private.base_set._get_version()
        
        base_set = self._private.base_set
        overlay_set = self._private.overlay_set
        
        base_set_keys = base_set.key
        overlay_set_keys = overlay_set.key
        removed_indices = []
        index = 0
        index_in_overlay = 0
        while index < len(base_set_keys) and index_in_overlay < len(overlay_set_keys):
            key1 = base_set_keys[index]
            key2 = overlay_set_keys[index_in_overlay]
            if key2 == key1:
                index_in_overlay += 1
                index += 1
            else:
                while key2 != key1:
                    removed_indices.append(index_in_overlay)
                    index_in_overlay += 1
                    if index_in_overlay >= len(overlay_set_keys):
                        break
                    key2 = overlay_set_keys[index_in_overlay]
            
        added_keys = []
        if index_in_overlay >= len(overlay_set_keys):
            added_keys = base_set_keys[index:]
        
        overlay_set.remove_particles_from_store(removed_indices)
        if len(added_keys) > 0:
            overlay_set.add_particles_to_store(added_keys)

    def can_extend_attributes(self):
        return self._private.overlay_set.can_extend_attributes()

    def __len__(self):
        return len(self._private.overlay_set)

    def _get_version(self):
        return self._private.overlay_set._get_version() +  self._private.base_set._get_version()


    def __getitem__(self, index):
        self._ensure_updated_set_properties()
        
        keys = self.get_all_keys_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Particle(keys, self)

    def _split_attributes(self, attributes):
        inbase = set(self._private.base_set.get_attribute_names_defined_in_store())
        attributes_inbase = []
        attributes_inoverlay = []
        indices_inbase = []
        indices_inoverlay = []
        for i,x in enumerate(attributes):
            if x in inbase:
                attributes_inbase.append(x)
                indices_inbase.append(i)
            else:
                attributes_inoverlay.append(x)
                indices_inoverlay.append(i)
        return (attributes_inbase, indices_inbase), (attributes_inoverlay, indices_inoverlay)

    def _split_attributes_and_values(self, attributes, values):
        inbase = set(self._private.base_set.get_attribute_names_defined_in_store())
        attributes_inbase = []
        attributes_inoverlay = []
        values_inbase = []
        values_inoverlay = []
        for x,y in zip(attributes, values):
            if x in inbase:
                attributes_inbase.append(x)
                values_inbase.append(y)
            else:
                attributes_inoverlay.append(x)
                values_inoverlay.append(y)
        return (attributes_inbase, values_inbase), (attributes_inoverlay, values_inoverlay)

    def add_particles_to_store(self, keys, attributes = [], values = []):
        self._ensure_updated_set_properties()
        
        (
            (attributes_inbase, values_inbase),
            (attributes_inoverlay, values_inoverlay)
        ) = self._split_attributes_and_values(attributes, values)


        self._private.base_set.add_particles_to_store(keys, attributes_inbase, values_inbase)
        self._private.overlay_set.add_particles_to_store(keys, attributes_inoverlay, values_inoverlay)

        self._private.base_version = self._private.base_set._get_version()
        
    def remove_particles_from_store(self, indices):
        self._ensure_updated_set_properties()
        
        indices = numpy.asarray(indices)

        self._private.base_set.remove_particles_from_store(indices[...,0])
        self._private.overlay_set.remove_particles_from_store(indices[...,1])

        self._private.base_version = self._private.base_set._get_version()
        


    def get_values_in_store(self, indices, attributes):
        self._ensure_updated_set_properties()
        
        (
            (attributes_inbase, indices_inbase),
            (attributes_inoverlay, indices_inoverlay)
        ) = self._split_attributes(attributes)


        if indices is None:
            indices0 = self._private.base_set.get_all_indices_in_store()
            indices1 = self._private.overlay_set.get_all_indices_in_store()
        else:
            indices0 = []
            indices1 = []
            for i0, i1 in indices:
                indices0.append(i0)
                indices1.append(i1)
            indices0 = numpy.asarray(indices0, dtype='int64')
            indices1 = numpy.asarray(indices1, dtype='int64')

        result = [None] * len(attributes)
        if len(attributes_inbase) > 0:
            values_inbase = self._private.base_set.get_values_in_store(indices0, attributes_inbase)
            for i, value in zip(indices_inbase, values_inbase):
                result[i] = value

        if len(attributes_inoverlay) > 0:
            values_inoverlay = self._private.overlay_set.get_values_in_store(indices1, attributes_inoverlay)
            for i, value in zip(indices_inoverlay, values_inoverlay):
                result[i] = value

        return result

    def set_values_in_store(self, indices, attributes, values):
        self._ensure_updated_set_properties()
        
        (
            (attributes_inbase, values_inbase),
            (attributes_inoverlay, values_inoverlay)
        ) = self._split_attributes_and_values(attributes, values)


        if indices is None:
            indices0 = self._private.base_set.get_all_indices_in_store()
            indices1 = self._private.overlay_set.get_all_indices_in_store()
        else:
            indices0 = []
            indices1 = []
            for i0, i1 in indices:
                indices0.append(i0)
                indices1.append(i1)
            indices0 = numpy.asarray(indices0, dtype='int64')
            indices1 = numpy.asarray(indices1, dtype='int64')

        if len(attributes_inbase) > 0:
            self._private.base_set.set_values_in_store(indices0, attributes_inbase, values_inbase)
        if len(attributes_inoverlay) > 0:
            self._private.overlay_set.set_values_in_store(indices1, attributes_inoverlay, values_inoverlay)


    def get_attribute_names_defined_in_store(self):
        result = list(self._private.base_set.get_attribute_names_defined_in_store())
        result.extend(self._private.overlay_set.get_attribute_names_defined_in_store())
        return result

    def get_all_keys_in_store(self):
        self._ensure_updated_set_properties()
        
        return self._private.overlay_set.get_all_keys_in_store()

    def get_all_indices_in_store(self):
        self._ensure_updated_set_properties()
        
        indices0 = self._private.base_set.get_all_indices_in_store()
        indices1 = self._private.overlay_set.get_all_indices_in_store()

        return zip(indices0, indices1)

    def get_indices_of_keys(self, keys):
        self._ensure_updated_set_properties()
        
        indices0 = self._private.base_set.get_indices_of_keys(keys)
        indices1 = self._private.overlay_set.get_indices_of_keys(keys)

        return zip(indices0, indices1)

    def has_key_in_store(self, key):
        self._ensure_updated_set_properties()
        
        return self._private.overlay_set.has_key_in_store(key)

    def _original_set(self):
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

    def compressed(self):
        return ParticlesWithUnitsConverted(self._private.particles.compressed(), self._private.converter)

    def get_valid_particles_mask(self):
        return self._private.particles.get_valid_particles_mask()

    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]

        if keys is ma.masked:
            return None
        elif hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Particle(keys, self)

    def _get_version(self):
        return self._private.particles._get_version()

    def shallow_copy(self):
        copiedParticles =  self._private.particles.shallow_copy()
        return ParticlesWithUnitsConverted(copiedParticles, self._private.converter)

    def unconverted_set(self):
        return self._private.particles


    def can_extend_attributes(self):
        return self._private.particles.can_extend_attributes()

    def add_particles_to_store(self, keys, attributes = [], values = []):
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_source_to_target(quantity)
            converted_values.append(converted_quantity)
        self._private.particles.add_particles_to_store(keys, attributes, converted_values)

    def remove_particles_from_store(self, keys):
        self._private.particles.remove_particles_from_store(keys)

    def get_values_in_store(self, indices, attributes):
        values = self._private.particles.get_values_in_store(indices, attributes)
        converted_values = []
        for quantity in values:
            if isinstance(quantity, LinkedArray):
                objects = quantity
                convert_objects = []
                for x in objects:
                    if x is None:
                        convert_objects.append(x)
                    else:
                        if isinstance(x, Particle):
                            convert_objects.append(ParticlesWithUnitsConverted(x.as_set(), self._private.converter)[0])
                        else:
                            convert_objects.append(ParticlesWithUnitsConverted(x, self._private.converter))
                convert_objects = LinkedArray(convert_objects)
                converted_values.append(convert_objects)
            else:
                converted_quantity = self._private.converter.from_target_to_source(quantity)
                converted_values.append(converted_quantity)
        return converted_values

    def set_values_in_store(self, indices, attributes, values):
        converted_values = []
        for quantity in values:
            converted_quantity = self._private.converter.from_source_to_target(quantity)
            converted_values.append(converted_quantity)
        self._private.particles.set_values_in_store(indices, attributes, converted_values)

    def get_attribute_names_defined_in_store(self):
        return self._private.particles.get_attribute_names_defined_in_store()

    def get_all_keys_in_store(self):
        return self._private.particles.get_all_keys_in_store()

    def get_all_indices_in_store(self):
        return self._private.particles.get_all_indices_in_store()

    def get_indices_of_keys(self, keys):
        return self._private.particles.get_indices_of_keys(keys)

    def has_key_in_store(self, key):
        return self._private.particles.has_key_in_store(key)

    def as_set(self):
        return ParticlesSubset(self, self.get_all_keys_in_store())

    def get_timestamp(self):
        timestamp = self._private.particles.get_timestamp()
        if not timestamp is None:
            timestamp = self._private.converter.from_target_to_source(timestamp)
        return timestamp

    def savepointsavepoint(self, timestamp=None):
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

    def get_unit_converter(self):
        return self._private.converter.generic_to_si



class ParticleInformationChannel(object):

    def __init__(self, from_particles, to_particles, attributes=None, target_names=None):
        self.from_particles = from_particles
        self.to_particles = to_particles
        self.attributes = attributes
        self.target_names = target_names
        self.from_version = -1
        self.to_version = -1
        self._reindex()

    def _reindex(self):
        if (
            (self.from_version == self.from_particles._get_version()) and
            (self.to_version == self.to_particles._get_version())
           ):
           return


        self.keys = self.intersecting_keys()
        self.from_indices = self.from_particles.get_indices_of_keys(self.keys)
        self.to_indices = self.to_particles.get_indices_of_keys(self.keys)
        self.from_version = self.from_particles._get_version()
        self.to_version = self.to_particles._get_version()

    def reverse(self):
        if self.target_names is None:
            attributes = self.attributes
            target_names = self.target_names
        else:
            attributes = self.target_names
            target_names = self.attributes

        return ParticleInformationChannel(
            self.to_particles,
            self.from_particles,
            attributes,
            target_names
        )

    def intersecting_keys(self):
        from_keys = self.from_particles.get_all_keys_in_store()
        to_keys = self.to_particles.get_all_keys_in_store()
        return numpy.intersect1d(from_keys,to_keys) #filter(lambda x : self.to_particles._has_key(x), from_keys)

    def copy_attributes(self, attributes, target_names = None):
        if target_names is None:
            target_names = attributes

        self._reindex()

        values = self.from_particles.get_values_in_store(self.from_indices, attributes)
        converted = []
        for x in values:
            if isinstance(x, LinkedArray):
                converted.append(x.copy_with_link_transfer(self.from_particles, self.to_particles))
            else:
                converted.append(x)

        self.to_particles.set_values_in_store(self.to_indices, target_names, converted)

    def copy(self):
        if not self.attributes is None:
            self.copy_attributes(self.attributes, self.target_names)
        elif not self.to_particles.can_extend_attributes():
            self.copy_overlapping_attributes()
        else:
            self.copy_all_attributes()

    def copy_all_attributes(self):
        names_to_copy = self.from_particles.get_attribute_names_defined_in_store()
        self.copy_attributes(list(names_to_copy))

    def copy_overlapping_attributes(self):
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
        data = self.from_particles.get_values_in_store(self.from_indices, [name,])
        self.to_particles.set_values_in_store(self.to_indices, [target_name,], data)


class ParticlesWithNamespacedAttributesView(AbstractParticleSet):
    """
    A view on prefixed attributes of a particle set.
    """

    def __init__(self, particles, namespace):
        AbstractParticleSet.__init__(self, particles)

        self._private.particles = particles
        self._private.namespace = namespace

    def compressed(self):
        return ParticlesWithNamespacedAttributesView(
            self._private.particles.compressed(),
            self._private.namespace
        )

    def get_valid_particles_mask(self):
        return self._private.particles.get_valid_particles_mask()

    def __getitem__(self, index):

        keys = self.get_all_keys_in_store()[index]

        if keys is ma.masked:
            return None
        elif hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Particle(keys, self)

    def _get_version(self):
        return self._private.particles._get_version()

    def shallow_copy(self):
        copiedParticles =  self._private.particles.shallow_copy()
        return ParticlesWithNamespacedAttributesView(
            copiedParticles,
            self._private.namespace
        )

    def unconverted_set(self):
        return self._private.particles

    def can_extend_attributes(self):
        return self._private.particles.can_extend_attributes()

    def add_particles_to_store(self, keys, attributes = [], values = []):
        namespace = self._private.namespace
        namespaced_attributes = [namespace + '__' + x for x in attributes]
        self._private.particles.add_particles_to_store(keys, namespaced_attributes, values)

    def remove_particles_from_store(self, keys):
        self._private.particles.remove_particles_from_store(keys)

    def get_values_in_store(self, indices, attributes):
        namespace = self._private.namespace
        namespaced_attributes = [namespace + '__' + x for x in attributes]
        return self._private.particles.get_values_in_store(indices, namespaced_attributes)


    def set_values_in_store(self, indices, attributes, values):
        namespace = self._private.namespace
        namespaced_attributes = [namespace + '__' + x for x in attributes]
        self._private.particles.set_values_in_store(indices, namespaced_attributes, values)

    def get_attribute_names_defined_in_store(self):
        names = self._private.particles.get_attribute_names_defined_in_store()
        namespace_prefix = self._private.namespace  + '__'
        len_namespace_prefix = len(namespace_prefix)
        return [x[len_namespace_prefix:] for x in names if x.startswith(namespace_prefix)]

    def get_all_keys_in_store(self):
        return self._private.particles.get_all_keys_in_store()

    def get_all_indices_in_store(self):
        return self._private.particles.get_all_indices_in_store()

    def get_indices_of_keys(self, keys):
        return self._private.particles.get_indices_of_keys(keys)

    def has_key_in_store(self, key):
        return self._private.particles.has_key_in_store(key)

    def as_set(self):
        return ParticlesSubset(self, self.get_all_keys_in_store())

    def get_timestamp(self):
        return self._private.particles.get_timestamp()

    def savepoint(self, timestamp=None):
        return ParticlesWithNamespacedAttributesView(
            self._private.particles.savepoint(timestamp),
            self._private.namespace
        )

    def previous_state(self):
        return ParticlesWithNamespacedAttributesView(
            self._private.particles.previous_state(),
            self._private.namespace
        )



class DomainAttribute(DerivedAttribute):
    """
    Combine multiple attributes into the same namespace
    """
    def  __init__(self, name):
        self.name = name

    def get_values_for_entities(self, instance):
        return ParticlesWithNamespacedAttributesView(instance, self.name)

    def set_values_for_entities(self, instance, value):
        raise AttributeError('"{0}" is already defined as a namespace attribute, you cannot assign a value to it'.format(self.name))

    def get_value_for_entity(self, instance, particle, index):
        namespaced_set = ParticlesWithNamespacedAttributesView(particle.particles_set, self.name)
        # or:
        # return namespaced_set[index]
        return Particle(
            particle.key,
            namespaced_set,
            particle._set_index,
            particle._set_version
        )

    def set_value_for_entity(self, instance, key, vector):
        raise AttributeError('"{0}" is already defined as a namespace attribute, you cannot assign a value to it'.format(self.name))




class Stars(Particles):
    pass

class Particle(object):
    """A physical object or a physical region simulated as a
    physical object (cloud particle).

    All attributes defined on a particle are specific for
    that particle (for example mass or position). A particle contains
    a set of attributes, some attributes are *generic* and applicable
    for multiple modules. Other attributes are *specific* and are
    only applicable for a single module.


    """
    __slots__ = ("key", "particles_set", "_set_index", "_set_version")

    # these are defined so that numpy conversion is way faster
    # otherwhise it would go through the __getattr__ function
    # which will slow it down by a factor 3
    if compare_version_strings(numpy.__version__, '1.7.0') < 0:
        __array_interface__ = {'shape':() }
    else:
        __array_interface__ = {'shape':(),'typestr':'|O4' }

    def __len__(self):
        raise AttributeError()
    def __iter__(self):
        raise AttributeError()

    __array_struct__ = UndefinedAttribute()
    __array__ = UndefinedAttribute()

    def __init__(self, key = None, particles_set = None, set_index = -1, set_version = -1, **keyword_arguments):
        if particles_set is None:
            if key == None:
                particles_set = Particles(1)
                key = particles_set.get_all_keys_in_store()[0]
            else:
                particles_set = Particles(1, keys = [key])

        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        object.__setattr__(self, "_set_index", set_index)
        object.__setattr__(self, "_set_version", set_version)

        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)

    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if self._set_index < 0 or self._set_version != self.particles_set._get_version():
            object.__setattr__(self, "_set_index", self.particles_set.get_indices_of_keys([self.key])[0])
            object.__setattr__(self, "_set_version", self.particles_set._get_version())

        self.particles_set._set_value_of_attribute(
            self._set_index,
            name_of_the_attribute,
            new_value_for_the_attribute
        )


    def __getattr__(self, name_of_the_attribute):
        try:
            if self._set_index < 0 or self._set_version != self.particles_set._get_version():
                object.__setattr__(self, "_set_index", self.particles_set.get_indices_of_keys([self.key])[0])
                object.__setattr__(self, "_set_version", self.particles_set._get_version())
            return self.particles_set._get_value_of_attribute(self, self._set_index, name_of_the_attribute)
        except Exception as ex:
            raise AttributeError("You tried to access attribute '{0}' but this attribute is not defined for this set.".format(name_of_the_attribute, ex))

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

    
    def copy(self):
        return self.particles_set.copy()._get_particle(self.key)
        
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
        if self._set_index < 0 or self._set_version != self.particles_set._get_version():
            object.__setattr__(self, "_set_index", self.particles_set.get_indices_of_keys([self.key])[0])
            object.__setattr__(self, "_set_version", self.particles_set._get_version())

        output = 'Particle('
        output += str(self.key)
        for name, value in self.particles_set._values_of_particle(self._set_index):
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

    def get_containing_set(self):
        return self.particles_set._original_set()

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




