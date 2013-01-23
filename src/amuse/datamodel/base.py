from amuse.support.core import CompositeDictionary
from amuse.support.core import OrderedDictionary
from amuse.support.core import late
from amuse.support import exceptions
from amuse.units import constants
from amuse.units import units
from amuse.units import quantities
from amuse.units.quantities import Quantity
from amuse.units.quantities import new_quantity
from amuse.units.quantities import is_quantity
from amuse.units.quantities import zero
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.units.quantities import as_vector_quantity

import numpy
import numpy.random
import random
import inspect
import warnings

class KeyGenerator(object):
    
    def next(self):
        pass
        
    def next_set_of_keys(self, length):
        pass
        
class BasicUniqueKeyGenerator(KeyGenerator):
    
    def __init__(self, lowest_unique_key = 1):
        self.lowest_unique_key = lowest_unique_key
    
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
        

class RandomNumberUniqueKeyGenerator(KeyGenerator):
    DEFAULT_NUMBER_OF_BITS = 64
    
    def __init__(self, number_of_bits = None):
        if number_of_bits is None:
            number_of_bits = self.DEFAULT_NUMBER_OF_BITS
        if number_of_bits > 64:
            raise exceptions.AmuseException("number of bits is larger than 64, this is currently unsupported!")
        self.number_of_bits = number_of_bits
        
        self.random = numpy.random.mtrand.RandomState()
        
    def next(self):
        return numpy.array([random.getrandbits(self.number_of_bits)], dtype=numpy.uint64)[0]
        
    def next_set_of_keys(self, length):
        if length == 0:
            return  []
        try:
            minint = -2** ((self.number_of_bits // 2) - 1)
            maxint = 2** ((self.number_of_bits // 2) - 1)
            low = self.random.random_integers(minint,maxint,length)
            high = self.random.random_integers(minint,maxint,length)
            return numpy.array(low + (high << 32), dtype=numpy.uint64)
        except:
             return numpy.array([random.getrandbits(self.number_of_bits) for i in range(length)], dtype='uint64')
       
        
UniqueKeyGenerator = RandomNumberUniqueKeyGenerator()

def set_sequential_key_generator(start_number):
    global UniqueKeyGenerator
    UniqueKeyGenerator = BasicUniqueKeyGenerator(lowest_unique_key = start_number)
    
def set_random_key_generator(number_of_bits = RandomNumberUniqueKeyGenerator.DEFAULT_NUMBER_OF_BITS):
    global UniqueKeyGenerator
    UniqueKeyGenerator = RandomNumberUniqueKeyGenerator(number_of_bits = number_of_bits)

class AttributeStorage(object):
    """
    Abstract base class of particle storage models.
    """
    __version__ = -1
    
    def add_particles_to_store(self, keys, attributes = [], values = []):
        pass
        
    def remove_particles_from_store(self, keys):
        pass
        
    def get_values_in_store(self, indices, attributes):
        pass
        
    def set_values_in_store(self, indices, attributes, list_of_values_to_set):
        pass
        
    def get_attribute_names_defined_in_store(self):
        pass
    
    def has_key_in_store(self, key):
        return False
        
    def get_all_keys_in_store(self):
        return []
        
    def get_all_indices_in_store(self):
        return []
        
    def __len__(self):
        return 0
        
    def get_value_in_store(self, index, attribute):
        return self.get_values_in_store([index],[attribute])[0][0]
    
    
class DerivedAttribute(object):
    """
    Abstract base class for calculated properties and 
    methods on sets.
    """
    def get_values_for_entities(self, particles):
        return None
    
    def set_values_for_entities(self, particles, value):
        raise exceptions.AmuseException("cannot set value of a DerivedAttribute")

    def get_value_for_entity(self, particles, particle, index):
        return None

    def set_value_for_entity(self, particles, key, value):
        raise exceptions.AmuseException("cannot set value of a DerivedAttribute")

class AbstractAttributeValue(object):
    def __str__(self):
        return self._values.__str__()
    
    def __repr__(self):
        return self._values.__repr__()
        
    def __eq__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values == other._values
        else:
            return self._values == other
    
    def __len__(self):
        return len(self.indices)
    
    def __add__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values + other._values
        else:
            return self._values + other
        
    def __radd__(self, other):
        return  self._values + other
        
    def __sub__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values - other._values
        else:
            return self._values - other
        
    def __rsub__(self, other):
        return other - self._values
    
    def __mul__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values * other._values
        else:
            return self._values * other
        
    def __rmul__(self, other):
        return self._values * other
    
    #def __imul__(self, other):
    #    if isinstance(other, Vector):
    #        newvalues = self._values * other._values
    #    else:
    #        newvalues = self._values * other
    #    self._set_values_for_entities(newvalues)
    #    return newvalues
    
    def __pow__(self, other):
        return self._values ** other

    def __div__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values / other._values
        else:
            return self._values / other
    
    def __rdiv__(self, other):
        return other / self._values
    
    def is_quantity(self):
        return hasattr(self._values, "is_quantity") and self._values.is_quantity()
    
    def __getattr__(self, name):
        return getattr(self._values, name)
    
     
    def __neg__(self):
        return -self._values
    
    def __abs__(self):
        return abs(self._values)
    
    def __lt__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values < other._values
        else:
            return self._values < other

    def __gt__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values > other._values
        else:
            return self._values > other
        
    def __ne__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values != other._values
        else:
            return self._values != other
        
    def __le__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values <= other._values
        else:
            return self._values <= other
    
    def __ge__(self, other):
        if isinstance(other, AbstractAttributeValue):
            return self._values >= other._values
        else:
            return self._values >= other
    
    def __dir__(self):
        return dir(self._values)
    
class VectorAttributeValue(AbstractAttributeValue):    
    def __init__(self, collection, indices, attribute_names):
        self.collection = collection
        self.indices = indices
        self.attribute_names = attribute_names
        self._values
            
    @late
    def _values(self):
        values = self.collection.get_values_in_store(self.indices, self.attribute_names)
          
        unit_of_the_values = None
        is_a_quantity = None
        for quantity in values:
            if unit_of_the_values is None:
                is_a_quantity = is_quantity(quantity)
                if is_a_quantity:
                    unit_of_the_values = quantity.unit
                break
    
        if is_a_quantity:
            results = []
            for quantity in values:
                if unit_of_the_values is None:
                    unit_of_the_values = quantity.unit
                results.append(quantity.value_in(unit_of_the_values))
        else:
            results = values
        
        results = numpy.array(results)
        for i in range(len(results.shape) - 1, 0, -1):
            results = numpy.swapaxes(results,0,i)
        
        if is_a_quantity:
            return unit_of_the_values.new_quantity(results)
        else:
            return results

    
    def __getitem__(self, index):
        
        return self._values.__getitem__(index)
    
    def __setitem__(self, index, value):
        raise Exception("must implement set item")
    
    def _set_values_for_entities(self, value):
        is_value_a_quantity = is_quantity(value)
        if is_value_a_quantity:
            vectors = value.number
        else:
            vectors = numpy.asanyarray(value)
            
        split = numpy.split(vectors, len(self.attribute_names), axis = vectors.ndim - 1)
        split = [x.reshape(x.shape[0] if len(x.shape) <= 2 else x.shape[:-1]) for x in split]
       
        if is_value_a_quantity:
            list_of_values = []
            for i in range(len(self.attribute_names)):
                values = value.unit.new_quantity(split[i])
                list_of_values.append(values)
        else:
            list_of_values = split
            
        self.collection.set_values_in_store(self.indices, self.attribute_names, list_of_values)
  
    
    
class VectorAttribute(DerivedAttribute):
    """
    Combine multiple attributes into a vecter attribute
    """
    def  __init__(self, attribute_names):
        self.attribute_names = attribute_names
    
    def get_values_for_entities(self, instance):
        #if 1:
        #    return VectorAttributeValue(instance, instance.get_all_indices_in_store(), self.attribute_names)
        
        values = instance.get_values_in_store(instance.get_all_indices_in_store(), self.attribute_names)
        unit_of_the_values = None
        is_a_quantity = None
        for quantity in values:
            if unit_of_the_values is None:
                is_a_quantity = is_quantity(quantity)
                if is_a_quantity:
                    unit_of_the_values = quantity.unit
                break
    
        if is_a_quantity:
            results = []
            for quantity in values:
                if unit_of_the_values is None:
                    unit_of_the_values = quantity.unit
                results.append(quantity.value_in(unit_of_the_values))
        else:
            results = values
        
        results = numpy.array(results)
        for i in range(len(results.shape) - 1, 0, -1):
            results = numpy.swapaxes(results,0,i)
        
        if is_a_quantity:
            return unit_of_the_values.new_quantity(results)
        else:
            return results

    def set_values_for_entities(self, instance, value):
        is_value_a_quantity = is_quantity(value)
        if is_value_a_quantity:
            vectors = value.number
        else:
            vectors = numpy.asanyarray(value)
            
        split = numpy.split(vectors, len(self.attribute_names), axis = vectors.ndim - 1)
        split = [x.reshape(x.shape[0] if len(x.shape) <= 2 else x.shape[:-1]) for x in split]
       
        if is_value_a_quantity:
            list_of_values = []
            for i in range(len(self.attribute_names)):
               values = value.unit.new_quantity(split[i])
               list_of_values.append(values)
        else:
            list_of_values = split
            
        instance.set_values_in_store(instance.get_all_indices_in_store(), self.attribute_names, list_of_values)
    
    def get_value_for_entity(self, instance, particle, index):
        values = instance._get_values_for_entity(index, self.attribute_names)
          
        unit_of_the_values = None
        is_a_quantity = None
        for quantity in values:
            if unit_of_the_values is None:
                is_a_quantity = is_quantity(quantity)
                if is_a_quantity:
                    unit_of_the_values = quantity.unit
                break
            
        if is_a_quantity:
            results = []
            for quantity in values:
                results.append(quantity.value_in(unit_of_the_values))
            return unit_of_the_values.new_quantity(results)
        else:
            return numpy.asarray(values)
            
    def set_value_for_entity(self, instance, key, vector):
        list_of_values = []
        for quantity in vector:
            if is_quantity(quantity):
                list_of_values.append(quantity.as_vector_with_length(1))
            else:
                list_of_values.append(quantity)
        instance._set_values_for_entity(key, self.attribute_names, list_of_values)


    
class CalculatedAttribute(DerivedAttribute):
    """
    Calculate the value of an attribute based
    on existing attributes.
    """
    def  __init__(self, function, attribute_names = None):
        self.function = function
        if attribute_names is None:
            arguments, varargs, kwargs, defaults = inspect.getargspec(function)
            self.attribute_names = arguments 
        else:
            self.attribute_names = attribute_names
    
    def get_values_for_entities(self, instance):
        values = instance.get_values_in_store(instance.get_all_indices_in_store(), self.attribute_names)
        return self.function(*values)
    
    def get_value_for_entity(self, particles,  particle, index):
        values = particles._get_values_for_entity(index, self.attribute_names)
        return self.function(*values)
        

class FunctionAttribute(DerivedAttribute):
    class BoundParticlesFunctionAttribute(object):
        def  __init__(self, function, particles):
            self.function = function
            self.particles = particles
            
        def __call__(self, *list_arguments, **keyword_arguments):
            return self.function(self.particles, *list_arguments, **keyword_arguments)
    
    class BoundParticleFunctionAttribute(object):
        def  __init__(self, function, particles, particle, index):
            self.function = function
            self.particles = particles
            self.index = index
            self.particle = particle
            
        def __call__(self, *list_arguments, **keyword_arguments):
            return self.function(self.particles, self.particle, *list_arguments, **keyword_arguments)
        
    
    def  __init__(self, particles_function = None, particle_function = None):
        self.particles_function = particles_function
        self.particle_function = particle_function
        
    def get_values_for_entities(self, particles):
        return self.BoundParticlesFunctionAttribute(self.particles_function, particles)
            
   
    def get_value_for_entity(self, particles, particle, index):
        return self.BoundParticleFunctionAttribute(self.particle_function, particles, particle, index)


class CollectionAttributes(object):
    """
    Objects of this class store attributes for
    a particles collection or a grid.
    
    These attributes are user set "meta" information such
    as a timestamp.
    """
    
    def __init__(self, attributes = None):
        if attributes is None:
            attributes = OrderedDictionary()
        else:
            attributes = attributes.copy()
        
        object.__setattr__(self, '_attributes', attributes)
    
    def __getattr__(self, name):
        try:
            return self._attributes[name]
        except KeyError:
            raise AttributeError('uknown attribute: {0!r}'.format(name))
        
    def __setattr__(self, name, value):
        self._attributes[name] = value
    
    def __getstate__(self):
        return self._attributes
        
    def __setstate__(self, data):
        object.__setattr__(self, '_attributes', data)
        
    def __str__(self):
        lines = []
        for name, value in self._attributes.iteritems():
            lines.append("{0}: {1}".format(name, value))
        return '\n'.join(lines)
    
    def _copy_for_collection(self, newcollection):
        return CollectionAttributes(self._attributes)
        

class PrivateProperties(object):
    """
    Defined for superclasses to store private properties.
    Every set has :meth:`__setattr__` defined.
    
    The :meth:`__setattr__` function will set all attributes
    of the entities in the set to the specified value(s).
    
    To be able to define attributes on the set itself we
    use an instance of this class, attributes can be 
    defined as::
    
        self._private.new_attribute = 'new value'
    
    Subclass implementers do not need to
    use the :meth:`object.__setattr__` syntax.
    
    For documentation about the :meth:`~object.__setattr__`
    call please see the 
    `python data model <http://docs.python.org/reference/datamodel.html>`_ 
    documentation on the python website.
    """
    pass

class AbstractSet(object):
    """
    Abstract superclass of all sets of particles and grids. 
    """
    GLOBAL_DERIVED_ATTRIBUTES = {}
        
    # this construct is needed to ensure that numpy
    # handles grids and particle sets as scalar objects
    # and not as sequences
    # if we put a grid in a numpy object array we want the
    # grid in a field of that array and not the contents of
    # the grid (i.e. the grid points)
    __array_interface__ = {'shape':()}
    
    #def __array_struct__(self):
    #    raise AttributeError()
    #def __array__(self):
    #    raise AttributeError()
    
    def __init__(self, original = None):
        if original is None:
            derived_attributes = self.GLOBAL_DERIVED_ATTRIBUTES
        else:
            derived_attributes = original._derived_attributes
            
        object.__setattr__(self, "_derived_attributes", CompositeDictionary(derived_attributes))
        object.__setattr__(self, "_private", PrivateProperties())
        
        self._private.collection_attributes = CollectionAttributes()
    
    @property
    def collection_attributes(self):
        return self._private.collection_attributes 
    
    @property
    def key(self):
        return self.get_all_keys_in_store()
        
    def __getattr__(self, name_of_the_attribute):
        if name_of_the_attribute == '__setstate__':
            raise AttributeError('type object {0!r} has no attribute {1!r}'.format(type(self), name_of_the_attribute))
        if name_of_the_attribute in self._derived_attributes:
            return self._derived_attributes[name_of_the_attribute].get_values_for_entities(self)
        else:
            try:
                return self._convert_to_entities_or_quantities(
                    self.get_all_values_of_attribute_in_store(name_of_the_attribute)
                )
            except Exception as ex:
                if name_of_the_attribute in self.get_attribute_names_defined_in_store():
                    raise
                else:
                    raise AttributeError("You tried to access attribute '{0}'"
                        " but this attribute is not defined for this set.".format(name_of_the_attribute))
    
    
    def check_attribute(self, value):
        if not (is_quantity(value) or hasattr(value, 'as_set')):
            try:
                return as_vector_quantity(value)
            except:
                return value
        else:
            return value
            
    def __setattr__(self, name_of_the_attribute, value):
        value = self.check_attribute(value)
        if name_of_the_attribute in self._derived_attributes:
            self._derived_attributes[name_of_the_attribute].set_values_for_entities(self, value)
        else:
            self.set_values_in_store(self.get_all_indices_in_store(), [name_of_the_attribute], [self._convert_from_entities_or_quantities(value)])
    
    def _get_value_of_attribute(self, particle, index, attribute):
        if attribute in self._derived_attributes:
            return self._derived_attributes[attribute].get_value_for_entity(self, particle, index)
        else:
            return self.get_value_in_store(index, attribute)
            
    def _get_values_for_entity(self, index, attributes):
        return [x[0] for x in self.get_values_in_store([index], attributes)]
        
    def _set_values_for_entity(self, index, attributes, values):
        return self.set_values_in_store([index], attributes, values)
    
    def _set_value_of_attribute(self, index, attribute, value):
        if attribute in self._derived_attributes:
            return self._derived_attributes[attribute].set_value_for_entity(self, index, value)
        else:
            return self.set_values_in_store(numpy.asarray([index]), [attribute], [value])
            
            
            
    def _convert_to_entities_or_quantities(self, x):
        if hasattr(x, 'unit') and x.unit.iskey():
            return self._subset(x.number)
        else:
            return x
        
    def _convert_from_entities_or_quantities(self, x):
        if isinstance(x, Quantity):
            return x 
        #elif isinstance(x, AbstractSet):
        #    return new_quantity( map(lambda y : (-1 if y is None else y.key), x), units.object_key) 
        else:
            return x
        
    #
    # Particle storage interface
    #
    
    def get_all_values_of_attribute_in_store(self, attribute):
        return self.get_values_in_store(self.get_all_indices_in_store(), [attribute])[0]
        
    def get_values_in_store(self, indices, attributes):
        pass
    
    def set_values_in_store(self, indices, attributes, values):
        pass
        
    def add_particles_to_store(self, indices, attributes, values):
        pass
        
    def remove_particles_from_store(self, indices):
        pass
    
    def get_all_keys_in_store(self):
        return []
        
    def get_all_indices_in_store(self):
        return []
        
    def has_key_in_store(self):
        return False
    
    def get_attribute_names_defined_in_store(self):
        return []

    def _original_set(self):
        return self
    
    def add_vector_attribute(self, name_of_the_attribute, name_of_the_components):
        self._derived_attributes[name_of_the_attribute] = VectorAttribute(name_of_the_components)
    
    @classmethod
    def add_global_vector_attribute(cls, name_of_the_attribute, name_of_the_components):
        """
        Define a *global* vector attribute, coupling two or more scalar attributes into
        one vector attribute. The vector will be defined for all particle sets
        created after calling this function.
        
        :argument name_of_the_attribute: Name to reference the vector attribute by. 
        :argument name_of_the_components: List of strings, each string a name of a scalar attribute.
        
        >>> from amuse.datamodel import Particles
        >>> particles = Particles(2)
        >>> Particles.add_global_vector_attribute('vel', ['vx','vy'])
        >>> particles.vx = [1.0 , 2.0] | units.m / units.s
        >>> particles.vy = [3.0 , 4.0] | units.m / units.s
        >>> particles.vel
        quantity<[[ 1.  3.], [ 2.  4.]] m / s>
        
        """
        cls.GLOBAL_DERIVED_ATTRIBUTES[name_of_the_attribute] = VectorAttribute(name_of_the_components)
    
    
    def add_calculated_attribute(self, name_of_the_attribute, function, attributes_names = None):
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
          
        
        >>> from amuse.datamodel import Particles
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
        
        self._derived_attributes[name_of_the_attribute] = CalculatedAttribute(function, attributes_names)
    
    
    @classmethod
    def add_global_calculated_attribute(cls, name_of_the_attribute, function, attributes_names = None):
        """
        Define a *global* vector attribute, coupling two or more scalar attributes into
        one vector attribute. The vector will be defined for all particle sets
        created after calling this function.
        
        :argument name_of_the_attribute: Name to reference the vector attribute by. 
        :argument function: Name of the function to call. 
        :argument attributes_names: List of strings, each string a name of a scalar attribute, if None uses argument names.
        
        
        >>> from amuse.datamodel import Particles
        >>> Particles.add_global_calculated_attribute("xy", lambda x, y : x * y)
        >>> particles = Particles(2)
        >>> particles.x = [1.0 , 2.0] | units.m 
        >>> particles.y = [3.0 , 4.0] | units.m
        >>> print particles.xy
        [3.0, 8.0] m**2
        >>> del Particles.GLOBAL_DERIVED_ATTRIBUTES['xy']
        
        """
        cls.GLOBAL_DERIVED_ATTRIBUTES[name_of_the_attribute] = CalculatedAttribute(function, attributes_names)
    
    
    def add_function_attribute(self, name_of_the_attribute, function, function_for_particle = None):
        """
        Define a function attribute, adding a function to the particles
        
        :argument name_of_the_attribute: Name to reference the vector attribute by. 
        :argument function: A function, first argument will be the particles.
        
        >>> from amuse.datamodel import Particles
        >>> particles = Particles(2)
        >>> particles.x = [1.0 , 2.0] | units.m
        >>> def sumx(p):
        ...   return p.x.sum()
        ...
        >>> particles.add_function_attribute("sum_of_x", sumx)
        >>> particles.sum_of_x()
        quantity<3.0 m>

        
        """
        
        self._derived_attributes[name_of_the_attribute] = FunctionAttribute(function, function_for_particle)
        
    
    @classmethod
    def add_global_function_attribute(cls, name_of_the_attribute, function, function_for_particle = None):
        """
        Define a function attribute, adding a function to the particles
        
        :argument name_of_the_attribute: Name to reference the attribute by. 
        :argument function: A function, first argument will be the particles.
        
        >>> from amuse.datamodel import Particles
        >>> def sumx(p):
        ...   return p.x.sum()
        ...
        >>> Particles.add_global_function_attribute("sum_of_x", sumx)
        >>> particles = Particles(2)
        >>> particles.x = [4.0 , 2.0] | units.m
        >>> particles.sum_of_x()
        quantity<6.0 m>
        >>> del Particles.GLOBAL_DERIVED_ATTRIBUTES['sum_of_x']
        
        """
        
        cls.GLOBAL_DERIVED_ATTRIBUTES[name_of_the_attribute] = FunctionAttribute(function, function_for_particle)
        
    
    def add_particle_function_attribute(self, name_of_the_attribute, function):
        """
        Define a function working on one particle
        
        :argument name_of_the_attribute: Name to reference the attribute by. 
        :argument function: A function, first argument will be the particle.
        
        >>> from amuse.datamodel import Particles
        >>> def xsquared(set, p):
        ...   return p.x * p.x
        ...
        >>> particles = Particles(2)
        >>> particles.add_particle_function_attribute("xsquared", xsquared)
        >>> particles.x = [4.0 , 2.0] | units.m
        >>> particles[0].xsquared()
        quantity<16.0 m**2>
        """
        self._derived_attributes[name_of_the_attribute] = FunctionAttribute(None, function)
        
    def __len__(self):
        return len(self.get_all_keys_in_store())

    def copy(self, memento = None, keep_structure = False):
        """ Creates a new in particle set and copies all attributes and 
        values into this set. 
        
        The history of the set is not copied over.
        
        
        Keyword arguments:
        memento -- internal, a dictionary to keep track of already copied sets in 
                   case of links between particles (default: None, will be created)
        keep_structure -- internal, if True, in case of a sub or super set make a copy of
                   the original set and return a new subset (default: False)
        """
        raise NotImplementedError()
        
    def copy_to_memory(self):
        warnings.warn("deprecated, pleasy change 'copy_to_memory' to 'copy'", DeprecationWarning)
        return self.copy()
    
    def _factory_for_new_collection(self):
        raise NotImplementedError()
        
    def copy_values_of_attribute_to(self, attribute_name, particles):
        """
        Copy values of one attribute from this set to the 
        other set. Will only copy values for the particles
        in both sets. See also :meth:`synchronize_to`.
        
        If you need to do this a lot, setup a dedicated
        channel.
        
        >>> from amuse.datamodel import Particles,Particle
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
        raise NotImplementedError()
        
    
    
    def __sub__(self, particles):
        """
        Returns a subset of the set without the given particle(s)
        Attribute values are not stored by the subset. The subset 
        provides a view on two or more sets of particles.
        
        :parameter particles: (set of) particle(s) to be subtracted from self.
        
        >>> from amuse.datamodel import Particles
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
        raise NotImplementedError()
    
    
    def add_particle(self, particle):
        """
        Add one particle to the set. 
        
        :parameter particle: particle to add
        
        >>> from amuse.datamodel import Particles,Particle
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
        
        >>> from amuse.datamodel import Particles
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
        
        >>> from amuse.datamodel import Particles
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
        
        >>> from amuse.datamodel import Particles, Particle
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
        
        >>> from amuse.datamodel import Particles
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
        
        >>> from amuse.datamodel import Particles
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
        
        >>> from amuse.datamodel import Particles
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
        values = map(lambda x: getattr(self, x), attributes)
        
        selections = selection_function(*values)
        selected_keys =  numpy.compress(selections, keys)
        
        return self._subset(selected_keys)
    
    def difference(self, other):
        """
        Returns a new subset containing the difference between
        this set and the provided set.
        
        >>> from amuse.datamodel import Particles
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
        
        >>> from amuse.datamodel import Particles,Particle
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
        raise NotImplementedError()
        
        
    def __dir__(self):
        """
        Utility function for introspection of paricle objects
        
        >>> from amuse.datamodel import Particles
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
        
    def stored_attributes(self):
        """
        Returns a list of the names of the attributes defined on
        the objects in this set (or grid).
        
        This list will not contain the set specific methods, or derived
        attributes (such as function attributes or vector attributes)
        for a list with all these attributes please use ``dir``.
        
        >>> from amuse.datamodel import Particles
        >>> particles = Particles(3)
        >>> particles.mass = [10.0, 20.0, 30.0] | units.kg
        >>> particles.x = [1.0, 2.0, 3.0] | units.m
        >>> print particles.stored_attributes()
        ['key', 'mass', 'x']
        """
        result = []
        result.append('key')
        result.extend(self.get_attribute_names_defined_in_store())
        return result
    
    def is_empty(self):
        return self.__len__()==0

    def get_value_in_store(self, key, attribute):
        return self.get_values_in_store(numpy.asarray([key]),[attribute])[0][0]
        
    def _convert_to_entity_or_quantity(self, x):
        if is_quantity(x):
            if x.unit.iskey():
                return self._subset([x.number])[0]
            else:
                return x
        else:
            return x

    def __add__(self, particles):
        """
        Returns a particle subset, composed of the given
        particle(s) and this particle set. Attribute values are
        not stored by the subset. The subset provides a view
        on two or more sets of particles.
        
        :parameter particles: (set of) particle(s) to be added to self.
        
        >>> from amuse.datamodel import Particles
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
        raise NotImplementedError()
    
    

    @classmethod
    def function_for_set(cls, function):
        cls.add_global_function_attribute(function.__name__, function)
        return function
    
    @classmethod
    def attribute_for_set(cls, function):
        cls.add_global_calculated_attribute(function.__name__, function)
        return function

    
    

class LinkedArray(numpy.ndarray):
    """Links between particles and particle sets are stored in LinkedArrays.
    """
    def __new__(cls, input_array):
        result = numpy.asarray(input_array).view(cls)
        return result

    def __array_finalize__(self, array_object):
        if array_object is None:
            return
        
    def copy(self, memento = None, keep_structure = False):
        from amuse.datamodel.particles import Particle
        from amuse.datamodel.grids import GridPoint
        
        if memento is None:
            memento = {}
            
        result = LinkedArray(self.flatten())
        index = 0
        for x in result:
            if x is None:
                result[index] = None
            elif isinstance(x, Particle):
                container = x.get_containing_set()
                if id(container) in memento:
                    copy_of_container = memento[id(container)]
                else:
                    copy_of_container = container.copy(memento, keep_structure)
                result[index] = copy_of_container._get_particle_unsave(x.key)
            elif isinstance(x, GridPoint):
                container = x.get_containing_set()
    
                if id(container) in memento:
                    copy_of_container = memento[id(container)]
                else:
                    copy_of_container = container.copy(memento, keep_structure)
                
                result[index] = GridPoint(x.index, copy_of_container)
            elif isinstance(x, AbstractSet):
                if id(x) in memento:
                    copy_of_container = memento[id(x)]
                else:
                    copy_of_container = x.copy(memento, keep_structure)
                result[index] = copy_of_container
            else:
                raise exceptions.AmuseException("unkown type in link {0}, copy not implemented".format(type(x)))
            index += 1
        
        return result.reshape(self.shape)
    
    def copy_with_link_transfer(self, from_container, to_container, must_copy = False, memento = None):
        from amuse.datamodel.particles import Particle
        from amuse.datamodel.grids import GridPoint
        
        if memento is None:
            memento = dict()
            
        result = LinkedArray(self.copy())
        index = 0
        for index in numpy.ndindex(*self.shape):
            if len(index) == 1:
                index = index[0]
            x = self[index]            
            if x is None:
                result[index] = None
            elif isinstance(x, Particle):
                container = x.get_containing_set()
                if from_container is None or container is from_container:
                    result[index] = to_container._get_particle_unsave(x.key)
                else:
                    result[index] = x
            elif isinstance(x, GridPoint):
                container = x.get_containing_set()
                if from_container is None or container is from_container:
                    result[index] = GridPoint(x.index, to_container)
                else:
                    result[index] = x
            elif isinstance(x, AbstractSet):
                if must_copy:
                    copy_of_container = x.copy(memento, keep_structure = True)
                    result[index] = copy_of_container
                else:
                    if from_container is None or x is from_container:
                        result[index] = to_container
                    else:
                        result[index] = x
            else:
                raise exceptions.AmuseException("unkown type in link {0}, transfer link not implemented".format(type(x)))
                
        return result
        
    def as_set(self):
        from amuse.datamodel.particles import Particle
        from amuse.datamodel.particles import ParticlesMaskedSubset
        
        linked_set = None
        mask = []
        keys = []
        index = 0
        
        flattened = LinkedArray(self.flatten())
        for x in flattened:
            if x is None:
                mask.append(True)
                keys.append(0)
            elif isinstance(x, Particle):
                original_set = x.as_set()._original_set()
                if linked_set is None:
                    linked_set = original_set
                elif not linked_set is original_set:
                    raise exceptions.AmuseException(
                        "could not convert the linked array to a subset as not all particles in the linked array are part of the same set"
                    )
                keys.append(x.key)
                mask.append(False)
            else:                        
                raise exceptions.AmuseException(
                    "could not convert the linked array to a subset as this array also contains sets of particles, grids or gridpoints"
                )
            index += 1
            
        if linked_set is None or len(linked_set) == 0:
           dtype = 'uint64'
        else:
           dtype = linked_set.get_all_keys_in_store().dtype
           
        masked_keys = numpy.ma.masked_array(
            numpy.asarray(
                keys,
                dtype = dtype
            ), 
            mask=mask
        )
        
        if linked_set is None:
            return ParticlesMaskedSubset(None, masked_keys)
        else:
            return linked_set._masked_subset(masked_keys)
    
    def to_print_list(self):
        from amuse.datamodel.particles import Particle
        
        result = []
        for x in self:
            if x is None:
                result.append('--')
            elif isinstance(x, Particle):
                result.append(x.key)
            else:
                result.append(type(x))
                
        return result
