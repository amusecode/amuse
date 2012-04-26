import numpy

from math import sqrt

from amuse.support import exceptions
from amuse.support.core import late
from amuse.support import console

from amuse.units.core import zero_unit

"""
"""
class Quantity(object):
    """
    A Quantity objects represents a scalar or vector with a
    specific unit. Quantity is an abstract base class
    for VectorQuantity and ScalarQuantity.

    Quantities should be constructed using *or-operator* ("|"),
    *new_quantity* or *unit.new_quantity*.

    Quantities emulate numeric types.

    Examples

    >>> from amuse.units import units
    >>> 100 | units.m
    quantity<100 m>
    >>> (100 | units.m) + (1 | units.km)
    quantity<1100.0 m>

    Quantities can be tested

    >>> from amuse.units import units
    >>> x = 100 | units.m
    >>> x.is_quantity()
    True
    >>> x.is_scalar()
    True
    >>> x.is_vector()
    False
    >>> v = [100, 200, 300] | units.g
    >>> v.is_quantity()
    True
    >>> v.is_scalar()
    False
    >>> v.is_vector()
    True

    Quantities can be converted to numbers

    >>> from amuse.units import units
    >>> x = 1000.0 | units.m
    >>> x.value_in(units.m)
    1000.0
    >>> x.value_in(units.km)
    1.0
    >>> x.value_in(units.g) # but only if the units are compatible!
    Traceback (most recent call last):
        File "<stdin>", line 1, in ?
    IncompatibleUnitsException: Cannot express m in g, the units do not have the same bases


    """
    __slots__ = ['unit']

    __array_priority__ = 101
    
    def __init__(self, unit):
        self.unit = unit

    def __str__(self):
        return console.current_printing_strategy.quantity_to_string(self)

    def is_quantity(self):
        """
        True for all quantities.
        """
        return True

    def is_scalar(self):
        """
        True for scalar quantities.
        """
        return False

    def is_vector(self):
        """
        True for vector quantities.
        """
        return False

    def __repr__(self):
        return 'quantity<'+str(self)+'>'

    def __add__(self, other):
        if self.unit.is_zero():
            return new_quantity(other.number, other.unit)
        else:
            other_in_my_units = to_quantity(other).as_quantity_in(self.unit)
            return new_quantity(self.number + other_in_my_units.number, self.unit)
    __radd__ = __add__

    def __sub__(self, other):
        if self.unit.is_zero():
            return -other
        else:
            other_in_my_units = to_quantity(other).as_quantity_in(self.unit)
            return new_quantity(self.number - other_in_my_units.number, self.unit)

    def __rsub__(self, other):
        if self.unit.is_zero():
            return new_quantity(other.number, other.unit)
        other_in_my_units = to_quantity(other).as_quantity_in(self.unit)
        return new_quantity(other_in_my_units.number - self.number, self.unit)

    def __mul__(self, other):
        other = to_quantity(other)
        return new_quantity_nonone(self.number * other.number, (self.unit * other.unit).to_simple_form())
        
    __rmul__ = __mul__

    def __pow__(self, other):
        return new_quantity(self.number ** other, self.unit ** other)

    def __div__(self, other):
        other = to_quantity(other)
        return new_quantity_nonone(self.number / other.number, (self.unit / other.unit).to_simple_form())

    def __rdiv__(self, other):
        return new_quantity_nonone(other / self.number, (1.0 / self.unit).to_simple_form())

    def in_(self, x):
        return self.as_quantity_in(x)

    def in_base(self):
        unit=self.unit.base_unit()
        return self.as_quantity_in(unit)

    def sqrt(self):
        """Calculate the square root of each component

        >>> from amuse.units import units
        >>> s1 = 144.0 | units.m**2
        >>> s1.sqrt()
        quantity<12.0 m>
        >>> v1 = [16.0, 25.0, 36.0] | units.kg
        >>> v1.sqrt()
        quantity<[4.0, 5.0, 6.0] kg**0.5>
        """
        return new_quantity(numpy.sqrt(self.number), (self.unit ** 0.5).to_simple_form())




    def as_quantity_in(self, another_unit):
        """
        Reproduce quantity in another unit.
        The new unit must have the same basic si quantities.

        :argument another_unit: unit to convert quantity to
        :returns: quantity converted to new unit
        """
        value_of_unit_in_another_unit = self.unit.as_quantity_in(another_unit)
        return new_quantity(self.number * value_of_unit_in_another_unit.number, another_unit)

    def value_in(self, unit):
        """
        Return a numeric value (for scalars) or array (for vectors)
        in the given unit.

        A number is returned without any unit information. Use this
        function only to transfer values to other libraries that have
        no support for quantities (for example plotting).

        :argument unit: wanted unit of the value
        :returns: number in the given unit

        >>> from amuse.units import units
        >>> x = 10 | units.km
        >>> x.value_in(units.m)
        10000.0

        """
        value_of_unit_in_another_unit = self.unit.value_in(unit)
        if value_of_unit_in_another_unit == 1.0:
            return self.number
        else:
            return self.number * value_of_unit_in_another_unit

    def __abs__(self):
        """
        Return the absolute value of this quantity

        >>> from amuse.units import units
        >>> x = -10 | units.km
        >>> print abs(x)
        10 km
        """
        return new_quantity(abs(self.number), self.unit)

    def __neg__(self):
        """
        Unary minus.

        >>> from amuse.units import units
        >>> x = -10 | units.km
        >>> print -x
        10 km
        """
        return new_quantity(-self.number, self.unit)

    def __lt__(self, other):
        return self.value_in(self.unit) < to_quantity(other).value_in(self.unit)

    def __gt__(self, other):
        return self.value_in(self.unit) > to_quantity(other).value_in(self.unit)

    def __eq__(self, other):
        return self.value_in(self.unit) == to_quantity(other).value_in(self.unit)

    def __ne__(self, other):
        return self.value_in(self.unit) != to_quantity(other).value_in(self.unit)

    def __le__(self, other):
        return self.value_in(self.unit) <= to_quantity(other).value_in(self.unit)

    def __ge__(self, other):
        return self.value_in(self.unit) >= to_quantity(other).value_in(self.unit)


class ScalarQuantity(Quantity):
    """
    A ScalarQuantity object represents a physical scalar
    quantity.
    """
    __slots__ = ['number']

    def __init__(self, number, unit):
        Quantity.__init__(self, unit)
        self.number = number

    def is_scalar(self):
        return True

    def as_vector_with_length(self, length):
        return VectorQuantity(numpy.ones(length, dtype=self.unit.dtype) * self.number, self.unit)


    def reshape(self, shape):
        if shape == -1 or (len(shape) == 1 and shape[0] == 1):
            return VectorQuantity([self.number], self.unit)
        else:
            raise exceptions.AmuseException("Cannot reshape a scalar to vector of shape '{0}'".format(shape))

    def __getitem__(self, index):
        if index == 0:
            return self
        else:
            raise exceptions.AmuseException("ScalarQuantity does not support indexing")

    def copy(self):
        return new_quantity(self.number, self.unit)

    def to_unit(self):
        in_base=self.in_base()
        return in_base.number * in_base.unit


    def __getstate__(self):
        return (self.unit, self.number)


    def new_zeros_array(self, length):
        array = numpy.zeros(length, dtype=self.unit.dtype)
        return new_quantity(array, self.unit)

    def __setstate__(self, tuple):
        self.unit = tuple[0]
        self.number = tuple[1]

    def sum(self, axis=None, dtype=None, out=None):
        return self
    def prod(self, axis=None, dtype=None):
        return self
    def amin(self, axis = None):
        return self
    def amax(self, axis = None):
        return self
    def sorted(self):
        return self
        
class VectorQuantity(Quantity):
    """
    A VectorQuantity object represents a physical vector
    quantity.

    >>> from amuse.units import units
    >>> v1 = [0.0, 1.0, 2.0] | units.kg
    >>> v2 = [2.0, 4.0, 6.0] | units.kg
    >>> v1 + v2
    quantity<[2.0, 5.0, 8.0] kg>
    >>> len(v1)
    3
    """
    __slots__ = ['_number']

    def __init__(self, array, unit):
        Quantity.__init__(self, unit)
        if unit is None:
            self._number = numpy.array((), dtype='float64')
        else:
            self._number = numpy.asarray(array, dtype=unit.dtype)

    @classmethod
    def new_from_scalar_quantities(cls, *values):
        unit = None
        array = []
        for x in values:
            if unit is None:
                unit = x.unit
            array.append(x.value_in(unit))
        return cls(array, unit)

    def aszeros(self):
        return new_quantity(numpy.zeros(self.shape, dtype=self.number.dtype), self.unit)

    
    def new_zeros_array(self, length):
        array = numpy.zeros(length, dtype=self.unit.dtype)
        return type(self)(array, self.unit)
        
    @classmethod
    def zeros(cls, length, unit):
        array = numpy.zeros(length, dtype=unit.dtype)
        return cls(array, unit)

    @classmethod
    def arange(cls, begin, end, step):
        unit = begin.unit
        #insert code here to throw exception if others have different unit
        
        val_begin = begin.value_in(unit)
        val_end = end.value_in(unit)
        val_step = step.value_in(unit)
        array = unit(numpy.arange(val_begin, val_end, val_step))
        return array

    @property
    def shape(self):
        return self.number.shape

    def flatten(self):
        return new_quantity(self.number.flatten(), self.unit)

    def is_vector(self):
        return True


    def as_vector_with_length(self, length):
        if length != len(self):
            raise exceptions.AmuseException("Can only return a vector with the same length")
        return self

    def as_vector_quantity(self):
        return self

    def __len__(self):
        return len(self._number)

    def sum(self, axis=None, dtype=None, out=None):
        """Calculate the sum of the vector components

        >>> from amuse.units import units
        >>> v1 = [0.0, 1.0, 2.0] | units.kg
        >>> v1.sum()
        quantity<3.0 kg>
        """
        return new_quantity(self.number.sum(axis, dtype, out), self.unit)

    def prod(self, axis=None, dtype=None):
        """Calculate the product of the vector components

        >>> from amuse.units import units
        >>> v1 = [1.0, 2.0, 3.0] | units.m
        >>> v1.prod()
        quantity<6.0 m**3>
        >>> v1 = [[2.0, 3.0], [2.0, 4.0], [5.0,3.0] ] | units.m
        >>> v1.prod()
        quantity<720.0 m**6>
        >>> v1.prod(0)
        quantity<[20.0, 36.0] m**3>
        >>> v1.prod(1)
        quantity<[6.0, 8.0, 15.0] m**2>

        >>> v1 = [[[2.0, 3.0], [2.0, 4.0]],[[5.0, 2.0], [3.0, 4.0]]] | units.m
        >>> v1.prod()
        quantity<5760.0 m**8>
        >>> v1.prod(0)
        quantity<[[ 10.   6.], [  6.  16.]] m**2>
        >>> v1.prod(1)
        quantity<[[  4.  12.], [ 15.   8.]] m**2>
        >>> v1.prod(2)
        quantity<[[ 6.  8.], [ 10.  12.]] m**2>
        """
        if axis == None:
            return new_quantity_nonone(self.number.prod(axis, dtype), self.unit ** numpy.prod(self.number.shape))
        else:
            return new_quantity_nonone(self.number.prod(axis, dtype), self.unit ** self.number.shape[axis])

    def length_squared(self):
        """Calculate the squared length of the vector.

        >>> from amuse.units import units
        >>> v1 = [2.0, 3.0, 4.0] | units.m
        >>> v1.length_squared()
        quantity<29.0 m**2>
        """
        return (self * self).sum()

    def length(self):
        """Calculate the length of the vector.

        >>> from amuse.units import units
        >>> v1 = [0.0, 3.0, 4.0] | units.m
        >>> v1.length()
        quantity<5.0 m>
        """
        return self.length_squared().sqrt()

    def lengths(self):
        """Calculate the length of the vectors in this vector.

        >>> from amuse.units import units
        >>> v1 = [[0.0, 3.0, 4.0],[2.0 , 2.0 , 1.0]] | units.m
        >>> v1.lengths()
        quantity<[5.0, 3.0] m>
        """
        return self.lengths_squared().sqrt()

    def lengths_squared(self):
        """Calculate the length of the vectors in this vector

        >>> from amuse.units import units
        >>> v1 = [[0.0, 3.0, 4.0],[4.0, 2.0, 1.0]] | units.m
        >>> v1.lengths_squared()
        quantity<[25.0, 21.0] m**2>
        """
        return (self.unit**2).new_quantity((self.number * self.number).sum(self.number.ndim - 1))

    def __getitem__(self, index):
        """Return the "index" component as a quantity.

        :argument index: index of the component, valid values
            for 3 dimensional vectors are: ``[0,1,2]``
        :returns: quantity with the same units

        >>> from amuse.units import si
        >>> vector = [0.0, 1.0, 2.0] | si.kg
        >>> print vector[1]
        1.0 kg
        >>> print vector[0:2]
        [0.0, 1.0] kg
        >>> print vector[[0,2,]]
        [0.0, 2.0] kg
        """

        return new_quantity( self._number[index], self.unit )

    def take(self, indices):
        return VectorQuantity(self._number.take(indices), self.unit)

    def put(self, indices, vector):
        try:
            if self.unit.is_zero():
                self.unit = vector.unit
            
            self._number.put(indices, vector.value_in(self.unit))
        except AttributeError:
            if not is_quantity(vector):
                raise ValueError("Tried to put a non quantity value in a quantity")
            else:
                raise

    def __setitem__(self, index, quantity):
        """Update the "index" component to the specified quantity.

        :argument index: index of the component, valid values
            for 3 dimensional vectors are: ``[0,1,2]``
        :quantity: quantity to set, will be converted to
            the unit of this vector

        >>> from amuse.units import si
        >>> vector = [0.0, 1.0, 2.0] | si.kg
        >>> g = si.kg / 1000
        >>> vector[1] = 3500 | g
        >>> print vector
        [0.0, 3.5, 2.0] kg
        """
        if self.unit.is_zero():
            self.unit = quantity.unit
        self._number[index] = quantity.value_in(self.unit)

    @property
    def number(self):
        return self._number

    @property
    def x(self):
        """The x axis component of a 3 dimensional vector.
        This is equavalent to the first component of vector.

        :returns: x axis component as a quantity

        >>> from amuse.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> print vector.x
        1.0 kg
        """
        return new_quantity(self.number[numpy.newaxis, ..., 0][0], self.unit)

    @property
    def y(self):
        """The y axis component of a 3 dimensional vector.
        This is equavalent to the second component of vector.

        :returns: y axis component as a quantity

        >>> from amuse.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> print vector.y
        2.0 kg
        """
        return new_quantity(self.number[numpy.newaxis, ..., 1][0], self.unit)

    @property
    def z(self):
        """The z axis component of a 3 dimensional vector.
        This is equavalent to the third component of vector.

        :returns: z axis component as a quantity

        >>> from amuse.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> print vector.z
        3.0 kg
        """
        return new_quantity(self.number[numpy.newaxis, ..., 2][0], self.unit)

    def indices(self):
        for x in len(self._number):
            yield x

    def copy(self):
        return new_quantity(self.number.copy(), self.unit)

    def norm_squared(self):
        return self.length_squared()

    def norm(self):
        return self.length()


    def append(self, scalar_quantity):
        """
        Append a scalar quantity to this vector.

        >>> from amuse.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> vector.append(4.0 | si.kg)
        >>> print vector
        [1.0, 2.0, 3.0, 4.0] kg
        """
        append_number = numpy.array(scalar_quantity.value_in(self.unit))
        # The following lines make sure that appending vectors works as expected,
        # e.g. ([]|units.m).append([1,2,3]|units.m) -> [[1,2,3]] | units.m
        # e.g. ([[1,2,3]]|units.m).append([4,5,6]|units.m) -> [[1,2,3],[4,5,6]] | units.m
        if (append_number.shape and (len(self._number) == 0 or self._number.shape[1:] == append_number.shape)):
            new_shape = [1 + self._number.shape[0]] + list(append_number.shape)
        else:
            new_shape = -1
        self._number = numpy.append(self._number, append_number).reshape(new_shape)

    def extend(self, vector_quantity):
        """
        Concatenate the vector quantity to this vector.
        If the units differ, the vector_quantity argument
        is converted to the units of this vector.

        >>> from amuse.units import units
        >>> vector1 = [1.0, 2.0, 3.0] | units.kg
        >>> vector2 = [1500, 2500, 6000] | units.g
        >>> vector1.extend(vector2)
        >>> print vector1
        [1.0, 2.0, 3.0, 1.5, 2.5, 6.0] kg
        """
        self._number = numpy.concatenate((self._number, vector_quantity.value_in(self.unit)))

    def prepend(self, scalar_quantity):
        """
        Prepend the scalar quantity before this vector.
        If the units differ, the scalar_quantity argument
        is converted to the units of this vector.

        >>> from amuse.units import units
        >>> vector1 = [1.0, 2.0, 3.0] | units.kg
        >>> vector1.prepend(0.0 | units.kg)
        >>> print vector1
        [0.0, 1.0, 2.0, 3.0] kg
        """
        self._number = numpy.concatenate(([scalar_quantity.value_in(self.unit)], self._number))

    def min(self, other):
        """
        Return the minimum of self and the argument.

        >>> from amuse.units import si
        >>> v1 = [1.0, 2.0, 3.0] | si.kg
        >>> v2 = [0.0, 3.0, 4.0] | si.kg
        >>> v1.min(v2)
        quantity<[0.0, 2.0, 3.0] kg>

        """
        other_in_my_units = other.as_quantity_in(self.unit)
        is_smaller_than = self.number < other_in_my_units.number
        values = numpy.where(is_smaller_than, self.number, other_in_my_units.number)
        return VectorQuantity(values, self.unit)

    def max(self, other):
        """
        Return the maximum of self and the argument.

        >>> from amuse.units import si
        >>> v1 = [1.0, 2.0, 3.0] | si.kg
        >>> v2 = [0.0, 3.0, 4.0] | si.kg
        >>> v1.max(v2)
        quantity<[1.0, 3.0, 4.0] kg>
        """
        other_in_my_units = other.as_quantity_in(self.unit)
        is_larger_than = self.number > other_in_my_units.number
        values = numpy.where(is_larger_than, self.number, other_in_my_units.number)
        return VectorQuantity(values, self.unit)

    def amax(self, axis = None):
        """
        Return the maximum along an axis.

        >>> from amuse.units import si
        >>> v1 = [1.0, 2.0, 3.0] | si.kg
        >>> v1.amax()
        quantity<3.0 kg>
        """

        return self.unit.new_quantity(numpy.amax(self.number, axis = axis))

    def amin(self, axis = None):
        """
        Return the minimum value along an axis.

        >>> from amuse.units import si
        >>> v1 = [1.0, 2.0, 3.0] | si.kg
        >>> v1.amin()
        quantity<1.0 kg>
        """

        return self.unit.new_quantity(numpy.amin(self.number, axis = axis))

    def sorted(self):
        """
        Return a new vector with all items sorted.

        >>> from amuse.units import si
        >>> v1 = [3.0, 1.0, 2.0] | si.kg
        >>> v1.sorted()
        quantity<[1.0, 2.0, 3.0] kg>
        """
        sorted_values = numpy.sort(self.number)
        return VectorQuantity(sorted_values, self.unit)
        
    def argsort(self, **options):
        """
        Returns the indices that would sort an array.

        >>> from amuse.units import si
        >>> v1 = [3.0, 1.0, 2.0] | si.kg
        >>> v1.argsort()
        array([1, 2, 0])
        """
        return numpy.argsort(self.number, **options)


    def sorted_with(self, *others):
        """
        Return a new vector with all items sorted. Perform
        all the same move operations on the other vectors.

        :argument: kind, the sort method for supported kinds see
            the numpy.sort documentation

        >>> from amuse.units import si
        >>> v1 = [3.0, 1.0, 2.0] | si.kg
        >>> v2 = [2.0, 3.0, 2.0] | si.m
        >>> v3 = [1.0, 4.0, 5.0] | si.s
        >>> list(v1.sorted_with(v2, v3))
        [quantity<[1.0, 2.0, 3.0] kg>, quantity<[3.0, 2.0, 2.0] m>, quantity<[4.0, 5.0, 1.0] s>]
        """
        indices = numpy.lexsort([self.number])
        vectors = []
        vectors.append(self)
        for x in others:
            vectors.append(x)

        for x in vectors:
            yield VectorQuantity(numpy.take(x.number, indices), x.unit)

    def accumulate(self):
        return VectorQuantity(numpy.add.accumulate(self.number), self.unit)

    def reshape(self, shape):
        return VectorQuantity(self.number.reshape(shape), self.unit)

    def transpose(self, axes=None):
        return VectorQuantity(self.number.transpose(axes), self.unit)

    def mean(self, axis=None, dtype=None, out=None):
        return new_quantity(self.number.mean(axis, dtype, out), self.unit)

    def median(self, **kwargs):
        return new_quantity(numpy.median(self.number, **kwargs), self.unit)

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        return new_quantity(self.number.std(axis, dtype, out, ddof), self.unit)

    def cross(self, other, axisa=-1, axisb=-1, axisc=-1, axis=None):
        """
        Return the cross product of this vector quantity with the supplied vector (quantity).
        """
        other = to_quantity(other)
        return new_quantity_nonone(
            numpy.cross(self.number, other.number, axisa=axisa, axisb=axisb, axisc=axisc, axis=axis), 
            (self.unit * other.unit).to_simple_form()
        )


    def tanh(self, out = None):
        if not self.unit.is_none():
            raise Exception("only none unit support for tanh")
        return numpy.tanh(self.number) | self.unit
        
    def sin(self):
        if not self.unit.is_none():
            raise Exception("only none unit support for sin")
        return new_quantity(numpy.sin(self.number), self.unit)
        
    def cos(self):
        if not self.unit.is_none():
            raise Exception("only none unit support for sin")
        return new_quantity(numpy.cos(self.number), self.unit)

    def __getstate__(self):
        return (self.unit, self.number)

    def __setstate__(self, tuple):
        self.unit = tuple[0]
        self._number = tuple[1]


class ZeroQuantity(Quantity):
    """
    A ZeroQuantity object represents zero in all units and
    can be used as the start for summing up purposes.

    >>> from amuse.units import si
    >>> x = zero
    >>> x += 2.0 | si.kg
    >>> x
    quantity<2.0 kg>

    """
    
    
    def __init__(self):
        Quantity.__init__(self, zero_unit())
        self.base = ()
        self.factor = 1
        self.number = 0.0
        self.dtype = 'float64'

    def is_scalar(self):
        """
        True for scalar quantities.
        """
        return False

    def is_vector(self):
        """
        True for vector quantities.
        """
        return False
    
    def is_non_numeric(self):
        return False
        
    def iskey(self):
        return False

    def __str__(self):
        return "zero"

    def __add__(self, other):
        return other

    def __sub__(self, other):
        return -other


    def __mul__(self, other):
        return self


    def __pow__(self, other):
        return self

    def __rmul__(self, other):
        return self


    def __div__(self, other):
        return self

    def __rdiv__(self, other):
        return other/self.number

    def in_base(self):
        return self


    def new_zeros_array(self, length):
        from amuse.units.core import zero_unit
        array = numpy.zeros(length, dtype=self.dtype)
        return new_quantity(array, zero_unit())

    def sqrt(self):
        return self

    def as_quantity_in(self, another_unit):
        return new_quantity(self.number, another_unit)

    def value_in(self, unit):
        return self.number

    def __abs__(self):
        return self

    def __neg__(self):
        return self

    def as_vector_with_length(self, length):
        return self.new_zeros_array(length)

    def __reduce__(self):
        return "zero"
        
    def __lt__(self, other):
        return 0 < to_quantity(other).value_in(other.unit)

    def __gt__(self, other):
        return 0 > to_quantity(other).value_in(other.unit)

    def __eq__(self, other):
        return 0 == to_quantity(other).value_in(other.unit)

    def __ne__(self, other):
        return 0 != to_quantity(other).value_in(other.unit)

    def __le__(self, other):
        return 0 <= to_quantity(other).value_in(other.unit)

    def __ge__(self, other):
        return 0 >= to_quantity(other).value_in(other.unit)



zero = ZeroQuantity()

class NonNumericQuantity(Quantity):
    """
    A Non Numeric Quantity object represents a quantity without
    a physical meaning.

    These Quantity objects cannot be used in
    numeric operations (like addition or multiplication). Also,
    conversion to another unit is not possible.

    Examples are string quantities or enumerated value quantities.

    >>> from amuse.units.core import enumeration_unit
    >>> my_unit = enumeration_unit(
    ...     "x",
    ...     "x",
    ...     [1,3,4],
    ...     ["first", "second", "third"])
    ...
    >>> 3 | my_unit
    quantity<3 - second>
    >>> (3 | my_unit).value_in(my_unit)
    3

    """

    def __init__(self, value, unit):
        Quantity.__init__(self, unit)
        self.value = value
        if not unit.is_valid_value(value):
            raise exceptions.AmuseException("<{0}> is not a valid value for {1!r}".format(value, unit))

    def as_quantity_in(self, another_unit):
        if not another_unit == self.unit:
            raise exceptions.AmuseException("Cannot convert non-numeric quantities in to another unit")

        return new_quantity(self.value, another_unit)

    def value_in(self, unit):
        if not unit == self.unit:
            raise exceptions.AmuseException("Cannot convert non-numeric quantities in to another unit")

        return self.value

    def __str__(self):
        return self.unit.value_to_string(self.value)

    def __repr__(self):
        return 'quantity<'+str(self.value)+ ' - ' +str(self)+'>'


    def as_vector_with_length(self, length):
        return VectorQuantity(numpy.array([self.value] * length), self.unit)

    def as_vector_quantity(self):
        return VectorQuantity([self.value], self.unit)
        
    def new_zeros_array(self, length):
        array = numpy.zeros(length, dtype=self.unit.dtype)
        return new_quantity(array, self.unit)

class AdaptingVectorQuantity(VectorQuantity):
    """
    Adapting vector quanity objects will adapt their units to the
    first object added to the vector
    """

    def __init__(self, value = [], unit = None):
        VectorQuantity.__init__(self, value, unit)
        del self._number
        self._number_list = list(value)
        if unit is None:
            self.append = self.append_start
            self.prepend = self.prepend_start
        else:
            self.append = self.append_normal
            self.prepend = self.prepend_normal
    
    def __getattr__(self, attribute):
        if attribute == "_number":
            if self.unit is None:
                return numpy.array(self._number_list)
            else:
                return numpy.array(self._number_list, dtype=self.unit.dtype)
        else:
            raise AttributeError
    
    def append_start(self, quantity):
        self.unit = quantity.unit
        self._number_list.append(quantity.value_in(self.unit))
        self.append = self.append_normal
        self.prepend = self.prepend_normal

    def append_normal(self, quantity):
        self._number_list.append(quantity.value_in(self.unit))
    
    prepend_start = append_start
    
    def prepend_normal(self, quantity):
        self._number_list.insert(0, quantity.value_in(self.unit))

    def extend(self, quantity):
        for x in quantity:
            self.append(x)

    def __setitem__(self, index, quantity):
        quantity_in_my_units = quantity.as_quantity_in(self.unit)
        self._number_list[index] = quantity_in_my_units.number
        self._remove_cached_number()

    def __getitem__(self, index):
        return new_quantity(self.number[index], self.unit)

    def __str__(self):
        if self.unit is None:
            return str(self._number_list)
        else:
            return console.current_printing_strategy.quantity_to_string(self)

def new_quantity(value, unit):
    """Create a new Quantity object.

    :argument value: numeric value of the quantity, can be
        a number or a sequence (list or ndarray)
    :argument unit: unit of the quantity
    :returns: new ScalarQuantity or VectorQuantity object
    """
    if isinstance(value, list):
        return VectorQuantity(value, unit)
    if isinstance(value, tuple):
        return VectorQuantity(value, unit)
    if isinstance(value, numpy.ndarray):
        return VectorQuantity(value, unit)
    if unit.is_non_numeric():
        return NonNumericQuantity(value, unit)
    return ScalarQuantity(value, unit)

def new_quantity_nonone(value, unit):
    """Create a new Quantity object.

    :argument value: numeric value of the quantity, can be
        a number or a sequence (list or ndarray)
    :argument unit: unit of the quantity
    :returns: new ScalarQuantity or VectorQuantity object
    """
    if not unit.base:
         return value * unit.factor
    if isinstance(value, list):
        return VectorQuantity(value, unit)
    if isinstance(value, tuple):
        return VectorQuantity(value, unit)
    if isinstance(value, numpy.ndarray):
        return VectorQuantity(value, unit)
    if unit.is_non_numeric():
        return NonNumericQuantity(value, unit)
    return ScalarQuantity(value, unit)
    
def is_quantity(input):
    return hasattr(input, "is_quantity")

def is_unit(input):
    if hasattr(input, "base"):
        return True
    else:
        return False

def as_vector_quantity(value):
    if not is_quantity(value):
        if hasattr(value, "__iter__"):
            result = AdaptingVectorQuantity()
            for subvalue in value:
                result.append(as_vector_quantity(subvalue))
            return result
        else:
            raise Exception("Cannot convert '{0!r}' to a vector quantity".format(value))
    return value
            
    
def to_quantity(input):
    if is_quantity(input):
        return input
    else:
        from amuse.units.units import none
        return new_quantity(input, none)
