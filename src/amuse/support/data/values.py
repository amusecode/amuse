"""
"""
import numpy

class Quantity(object):
    """
    A Quantity objects represents a scalar or vector with a
    specific unit. Quantity is an abstract base class
    for VectorQuantity and ScalarQuantity.
    
    Quantities should be constructed using *or-operator* ("|"),
    *new_quantity* or *unit.new_quantity*.
    
    Quantities emulate numeric types.
    
    Examples
    
    >>> from amuse.support.units import units
    >>> 100 | units.m
    quantity<100 m>
    >>> (100 | units.m) + (1 | units.km)
    quantity<1100.0 m>
    
    Quantities can be tested
       
    >>> from amuse.support.units import units
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
    
    >>> from amuse.support.units import units
    >>> x = 1000 | units.m
    >>> x.value_in(units.m)
    1000.0
    >>> x.value_in(units.km)
    1.0
    >>> x.value_in(units.g) # but only if the units are compatible!
    Traceback (most recent call last):
        File "<stdin>", line 1, in ?
    Exception: Cannot expres: g in m

    
    """
    __slots__ = ['unit']
    
    def __init__(self, unit):
        self.unit = unit
        
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
        other_in_my_units = other.as_quantity_in(self.unit)
        return new_quantity(self.number + other_in_my_units.number , self.unit)
        
    def __sub__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return new_quantity(self.number - other_in_my_units.number , self.unit)
        

    def __mul__(self, other):
        if hasattr(other, "unit"):
            return new_quantity(self.number * other.number , (self.unit * other.unit).to_simple_form())
        else:
            return new_quantity(self.number * other , self.unit.to_simple_form())
        
    
    def __pow__(self, other):
        return new_quantity(self.number ** other , self.unit ** other)
        
    def __rmul__(self, other):
        return self.__mul__(other)
        
      
    def __div__(self, other):
        if hasattr(other, "unit"):
            return new_quantity(self.number / other.number , (self.unit / other.unit).to_simple_form())
        else:
            return new_quantity(self.number / other , self.unit.to_simple_form())
    
    def __rdiv__(self, other):
        return self.__div__(other)
            
    def in_(self, x):
        return self.as_quantity_in(x)
    
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
        
        >>> from amuse.support.units import units
        >>> x = 10 | units.km
        >>> x.value_in(units.m)
        10000.0
        
        """
        value_of_unit_in_another_unit = self.unit.value_in(unit)
        return self.number * value_of_unit_in_another_unit

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
      
    def __str__(self):
        unit_str = str(self.unit)
        if unit_str:
            return str(self.number) + ' ' + unit_str
        else:
            return str(self.number)
                                
    def __lt__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number < other_in_my_units.number
        
    def __gt__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number > other_in_my_units.number
        
    def __eq__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number == other_in_my_units.number
        
    def __neq__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number != other_in_my_units.number

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)
            
    def copy(self):
        return new_quantity(self.number, self.unit)
            
class VectorQuantity(Quantity):
    """
    A VectorQuantity object represents a physical vector 
    quantity.
    
    >>> from amuse.support.units import units
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
        self._number = numpy.array(array)
        
    @classmethod
    def new_from_scalar_quantities(cls, *values):
        unit = None
        array = []
        for x in values:
            if unit is None:
                unit = x.unit
            array.append(x.value_in(unit))
        return cls(array, unit)
    
    def is_vector(self):
        return True
    
    def __len__(self):
        return len(self._number)
    
    def sum(self):
        """Calculate the sum of the vector components

        >>> from amuse.support.units import units
        >>> v1 = [0.0, 1.0, 2.0] | units.kg
        >>> v1.sum()
        quantity<3.0 kg>
        """
        return new_quantity(numpy.sum(self.number), self.unit)
        
    def sqrt(self):
        """Calculate the square root of each component

        >>> from amuse.support.units import units
        >>> v1 = [16.0, 25.0, 36.0] | units.kg
        >>> v1.sqrt()
        quantity<[4.0, 5.0, 6.0] kg**0.5>
        """
        return new_quantity(numpy.sqrt(self.number), (self.unit ** 0.5).to_simple_form())
    
    def length(self):
        """Calculate the length of the vector.
        
        >>> from amuse.support.units import units
        >>> v1 = [0.0, 3.0, 4.0] | units.m
        >>> v1.length()
        quantity<5.0 m>
        """
        squared = self * self
        dot_product = squared.sum()
        return (dot_product ** 0.5).as_quantity_in(self.unit)
        
    def __getitem__(self, index):
        """Return the "index" component as a quantity.
        
        :argument index: index of the component, valid values
            for 3 dimensional vectors are: ``[0,1,2]``
        :returns: quantity with the same units
        
        >>> from amuse.support.units import si
        >>> vector = [0.0, 1.0, 2.0] | si.kg
        >>> print vector[1]
        1.0 kg
        """
        
        return new_quantity( self._number[index] , self.unit )
        
    
    def __setitem__(self, index, quantity):
        """Update the "index" component to the specified quantity.
        
        :argument index: index of the component, valid values
            for 3 dimensional vectors are: ``[0,1,2]``
        :quantity: quantity to set, will be converted to
            the unit of this vector
        
        >>> from amuse.support.units import si
        >>> vector = [0.0, 1.0, 2.0] | si.kg
        >>> g = si.kg / 1000
        >>> vector[1] = 3500 | g
        >>> print vector
        [0.0, 3.5, 2.0] kg
        """
        quantity_in_my_units = quantity.as_quantity_in(self.unit)
        self._number[index] = quantity_in_my_units.number
        
    @property
    def number(self):
        return self._number
        
    @property
    def x(self):
        """The x axis component of a 3 dimensional vector.
        This is equavalent to the first component of vector.
        
        :returns: x axis component as a quantity
        
        >>> from amuse.support.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> print vector.x
        1.0 kg
        """
        return self[0]
        
    @property
    def y(self):
        """The y axis component of a 3 dimensional vector.
        This is equavalent to the second component of vector.
        
        :returns: y axis component as a quantity
        
        >>> from amuse.support.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> print vector.y
        2.0 kg
        """
        return self[1]
        
    @property
    def z(self):
        """The z axis component of a 3 dimensional vector.
        This is equavalent to the third component of vector.
        
        :returns: z axis component as a quantity
        
        >>> from amuse.support.units import si
        >>> vector = [1.0, 2.0, 3.0] | si.kg
        >>> print vector.z
        3.0 kg
        """
        return self[2]
        
    def __str__(self):
        unit_str = str(self.unit)
        array_str = '[' + ', '.join([str(x) for x in self._number]) + ']'
        if unit_str:
            return array_str + ' ' + unit_str
        else:
            return array_str
            

    def __lt__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number < other_in_my_units.number
        
    def __gt__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number > other_in_my_units.number
        
    def __eq__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        if not len(self.number) == len(other_in_my_units.number):
            return False
            
        return all(
            map(
                lambda x,y : x == y, 
                self.number, 
                other_in_my_units.number
                )
            )
        
    def __neq__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number != other_in_my_units.number   
        
    def indices(self):
        for x in len(self._number):
            yield x
            
    def copy(self):
        return new_quantity(self.number.copy(), self.unit)
                 

class NonNumericQuantity(Quantity):
    """
    A Non Numeric Quantity object represents a quantity without
    a physical meaning. 
    
    These Quantity objects cannot be used in 
    numeric operations (like addition or multiplication). Also,
    conversion to another unit is not possible.
    
    Examples are string quantities or enumerated value quantities.
    
    >>> from amuse.support.units.core import enumeration_unit
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
            raise Exception("<{0}> is not a valid value for {1!r}".format(value, unit))
    
    def as_quantity_in(self, another_unit): 
        if not another_unit == self.unit:
            raise Exception("Cannot convert non-numeric quantities in to another unit")
            
        return new_quantity(self.value, another_unit)

    def value_in(self, unit):
        if not unit == self.unit:
            raise Exception("Cannot convert non-numeric quantities in to another unit")
        
        return self.value
        
    def __str__(self):
        return self.unit.value_to_string(self.value)
        
    def __repr__(self):
        return 'quantity<'+str(self.value)+ ' - ' +str(self)+'>'
        
    def __eq__(self, other):
        if not other.unit == self.unit:
            return False
        return self.value == other.value
        
    
    
        
def new_quantity(value, unit):
    """Create a new Quantity object.
    
    :argument value: numeric value of the quantity, can be 
        a number or a sequence (list or ndarray)
    :argument unit: unit of the quantity
    :returns: new ScalarQuantity or VectorQuantity object
    """
    if isinstance(value, list):
        return VectorQuantity(value, unit)
    if isinstance(value, numpy.ndarray):
        return VectorQuantity(value, unit)
    if unit.is_non_numeric():
        return NonNumericQuantity(value, unit)
    return ScalarQuantity(value, unit)
    
    
