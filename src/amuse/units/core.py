from amuse.support.core import late
from amuse.support import exceptions
import numpy


from amuse.support.core import memoize

from amuse.support.core import MultitonMetaClass

class system(object):
    ALL = {}

    def __init__(self, name):
        self.name = name
        self.bases = []
        self.mapping_from_base_name_to_base = {}
        self.ALL[self.name] = self
        self.index = len(self.ALL)

    def reference_string(self):
        return "{0}.get({1!r})".format('system', self.name)

    def add_base(self, unit):
        unit.system = self
        unit.index = len(self.bases)
        self.bases.append(unit)
        self.mapping_from_base_name_to_base[unit.quantity] = unit

    def base(self, name):
        return self.mapping_from_base_name_to_base[name]

    @classmethod
    def get(cls, name):
        try:
            return cls.ALL[name]
        except KeyError as ex:
            from amuse.units import nbody_system
            from amuse.units import si
            return cls.ALL[name]


    def __reduce__(self):
        return (get_system_with_name, (self.name,))


    def __str__(self):
        return self.name


class unit(object):
    """
    Abstract base class for unit objects.
        
    Two classes of units are defined:

    base units
        The base units in a given system of units. For SI, these
        are meter, kilogram, second, ampere, kelvin, mole and 
        candele. See the si module :mod:`amuse.units.si`
    derived units
        Derived units are created by dividing or multiplying
        with a number or with another unit. For example, 
        to get a velocity unit we can devine vel = 1000 * m / s

    Units can also be named, by creating a named unit.
    """
    __array_priority__ = 100

    def __mul__(self, other):

        if isinstance(other, unit):
            return mul_unit(self, other)
        else:
            return other*self
#            return factor_unit(other, self)

    def __truediv__(self, other):
        if isinstance(other, unit):
            return div_unit(self, other)
        else:
            return (1.0/other)*self
#            return factor_unit(1.0 / other, self)

    def __rmul__(self, other):
        if other == 1:
            return self
        else:
            return factor_unit(other, self)

    def __ror__(self, value):
        """Create a new Quantity object.
    
        :argument value: numeric value of the quantity, can be 
            a number or a sequence (list or ndarray)
        :returns: new ScalarQuantity or VectorQuantity object 
            with this unit
            
        Examples
        
        >>> from amuse.units import units
        >>> 100 | units.kg
        quantity<100 kg>
        """
        return self.new_quantity(value)

    def __rtruediv__(self, other):
        return factor_unit(other, pow_unit(-1,self))

    def __div__(self, other):
        return self.__truediv__(other)

    def __rdiv__(self, other):
        return self.__rtruediv__(other)

    def __pow__(self, other):
        if other == 1:
            return self
        else:
            return pow_unit(other, self)

    def __call__(self, x):
        return self.new_quantity(x)

    def __eq__(self, other):
        if self is other:
            return True
        elif isinstance(other, unit):
            return self.base == other.base and self.factor == other.factor
        else:
            return False


    def __ne__(self, other):
        if isinstance(other, unit):
            if (isinstance(self, base_unit) and isinstance(other, base_unit)) or \
                    isinstance(self, nonnumeric_unit) or isinstance(other, nonnumeric_unit):
                return NotImplemented
            return self.base != other.base and self.factor != other.factor
        else:
            return True

    def __hash__(self):
        return self._hash

    @late
    def _hash(self):
        return hash(id(self))

    @property
    def dtype(self):
        return None

    @property
    def number(self):
        return 1.0

    @property
    def unit(self):
        return self

    def is_zero(self):
        return False

    def iskey(self):
        return False

    def new_quantity(self, value):
        """Create a new Quantity object.

        :argument value: numeric value of the quantity, can be 
            a number or a sequence (list or ndarray)
        :returns: new ScalarQuantity or VectorQuantity object 
            with this unit
        """
        from amuse.units import quantities
        return quantities.new_quantity(value, self)

    def to_simple_form(self):
        """Convert unit to a form with only one factor and powers
        
        :result: Unit with only a factor and power terms
        
        >>> from amuse.units import units
        >>> N = (units.m * units.kg) / (units.s * units.s)
        >>> N
        unit<m * kg / (s * s)>
        >>> J = N * units.m
        >>> J
        unit<m * kg / (s * s) * m>
        >>> J.to_simple_form()
        unit<m**2 * kg * s**-2>
        """

        if not self.base:
            return none_unit('none', 'none') * self.factor

        result = self.factor
        for n, base in self.base:
            if n == 1:
                if result == 1:
                    result = base
                else:
                    result = result * base
            else:
                result =  result * (base ** n)

        return result

    def to_reduced_form(self):
        """Convert unit to a reduced (simpler) form 
        """

        if not self.base:
            return none_unit('none', 'none') * self.factor

        total_factor = 1
        combined_unit = None
        for factor, power, unit in self.get_parts_with_power():
            total_factor *= factor
            if power == 0:
                pass
            else:
                if combined_unit is None:
                    combined_unit = unit ** power
                else:
                    combined_unit =  combined_unit * (unit ** power)
        if total_factor == 1:
            return combined_unit
        else:
            return factor_unit(total_factor , combined_unit)



    def to_factor_and_reduced_form(self):
        """Convert unit to a reduced (simpler) form 
        """

        if not self.base:
            return none_unit('none', 'none') * self.factor

        total_factor = 1
        combined_unit = None
        for factor, power, unit in self.get_parts_with_power():
            total_factor *= factor
            if power == 0:
                pass
            else:
                if combined_unit is None:
                    combined_unit = unit ** power
                else:
                    combined_unit =  combined_unit * (unit ** power)

        return total_factor , combined_unit

    def are_bases_equal(self, other):
        if len(self.base) != len(other.base):
            return False
        for n1, unit1 in sorted(self.base, key=lambda x: x[1].index):
            found = False
            for n2, unit2 in other.base:
                if unit1 == unit2:
                    if not n2 == n1:
                        return False
                    found = True
                    break
            if not found:
                return False
        return True

    def _compare_bases(self, other, eps = None):
        if len(self.base) != len(other.base):
            return False
        if eps is None:
            eps = numpy.finfo(numpy.double).eps
        for (n1, unit1), (n2, unit2) in zip(self.base, other.base):
            if not unit1 == unit2:
                return False
            if n1 == n2:
                continue
            else:
                if abs(n1 - n2) < eps:
                    continue
                if abs(n2) > abs(n1):
                    relativeError = abs((n1 - n2) * 1.0 / n2)
                else:
                    relativeError = abs((n1 - n2) * 1.0 / n1)
                if relativeError <= eps:
                    continue
                else:
                    return False
        return True

    @memoize
    def conversion_factor_from(self, x):
        if x.base is None:
            return self.factor * 1.0
        elif self._compare_bases(x):
            this_factor = self.factor * 1.0
            other_factor = x.factor
            return 1*(this_factor == other_factor) or this_factor / other_factor
        else:
            raise IncompatibleUnitsException(x, self)


    def in_(self, x):
        """Express this quantity in the given unit
        
        :argument unit: The unit to express this quantity in
        :result: A Quantity object
 
        Examples

        >>> from amuse.units import units
        >>> l = 1 | units.AU
        >>> l.in_(units.km)
        quantity<149597870.691 km>
        
        """

        return self.as_quantity_in(x)

    def as_quantity_in(self, unit):
        """Express this unit as a quantity in the given unit
        
        :argument unit: The unit to express this unit in
        :result: A Quantity object
        
        Examples
        
        >>> from amuse.units import units
        >>> ton = 1000 * units.kg
        >>> ton.as_quantity_in(units.kg)
        quantity<1000.0 kg>
        """
        from amuse.units import quantities
        if isinstance(unit, quantities.Quantity):
            raise exceptions.AmuseException("Cannot expres a unit in a quantity")

        factor = self.conversion_factor_from(unit)
        return quantities.new_quantity(factor, unit)

    def value_in(self, unit):
        """
        Return a numeric value of this unit in the given unit.
        Works only when the units are compatible, i.e. from
        tonnage to kg's.
        
        A number is returned without any unit information.
        
        :argument unit: wanted unit of the value
        :returns: number in the given unit
        
        >>> from amuse.units import units
        >>> x = units.km
        >>> x.value_in(units.m)
        1000.0
        
        """
        return self.conversion_factor_from(unit)

    def __repr__(self):
        return 'unit<'+str(self)+'>'

    def combine_bases(self, base1, base2):
        indexed1 = [None] * 7
        for n1, unit1 in base1:
            indexed1[unit1.index] = (n1, unit1)

        indexed2 = [None] * 7
        for n2, unit2 in base2:
            indexed2[unit2.index] = (n2, unit2)

        result = []
        for sub1, sub2 in zip(indexed1, indexed2):
            if not sub1 is None:
                if not sub2 is None:
                    if sub1[1] == sub2[1]:
                        result.append((sub1[0], sub2[0], sub1[1]))
                    else:
                        raise exceptions.AmuseException("Cannot combine units from "
                            "different systems: {0} and {1}".format(sub1[1], sub2[1]))
                else:
                    result.append((sub1[0], 0, sub1[1]))
            elif not sub2 is None:
                result.append((0, sub2[0], sub2[1]))
        return result

    def has_same_base_as(self, other):
        """Determine if the base of other is the same as the
        base of self.
        
        :argument other: unit to compare base to
        :result: True, if bases are compatiple.
        
        >>> from amuse.units import units
        >>> mps = units.m / units.s
        >>> kph = units.km / units.hour
        >>> mps.has_same_base_as(kph)
        True
        >>> mps.has_same_base_as(units.km)
        False
        
        """
        return other.base == self.base

    def base_unit(self):
        if not self.base:
            return none_unit('none', 'none')

        unit = 1
        for n, base in self.base:
            if n == 1:
                unit = unit*base 
            else:
                unit = unit*(base ** n)        
        return unit

    def is_non_numeric(self):
        return False

    def is_generic(self):
        return False

    def is_none(self):
        return False

    def get_parts_with_power(self):
        """
        The parts of this units as a list of tuple with factor, power and unit
        """
        return ((1.0, 1, self),)

    def convert_result_value(self, method, definition, value):
        return self.new_quantity(value)

    def convert_argument_value(self, method, definition, value):
        return value.value_in(self)

    def append_result_value(self, method, definition, value, result):
        result.append(self.convert_result_value(method, definition, value))


    def to_array_of_floats(self):
        """Represent a unit as an array of 8 64-bit floats. First float represents the factor, the other 7 the power of each base unit.
        Cannot be used for non numeric units
        """
        result = numpy.zeros(9, dtype=numpy.float64)
        if not self.base:
            return result

        result[0] = self.factor

        for n, base in self.base:
            result[base.index + 2] = n
            result[1] = base.system.index
        return result

    def describe_array_of_floats(self):
        """Create a human readable description of the array of floats
        """
        if not self.base:
            return 'not a numerical unit'

        parts = ['factor']
        parts.extend(['-']*8)

        for n, base in self.base:
            if n != 0:
                parts[base.index + 2] = str(base)
            else:
                parts[base.index + 2] = '-'

            parts[1] = str(base.system)
        return ', '.join(parts)



class base_unit(unit):
    """
    base_unit objects are  orthogonal, indivisable units 
    of a sytem of units.
    
    A system of units contains a set of base units
    
    :argument quantity: name of the base quantity, for example *length*
    :argument name: name of the unit, for example *meter*
    :argument symbol: symbol of the unit, for example *m*
    :argument system: system of units object
    
    >>> cgs = system("cgs")
    >>> cm = base_unit("length", "centimetre", "cm", cgs)
    >>> cm
    unit<cm>
    """
    def __init__(self, quantity, name, symbol, system):
        self.quantity = quantity
        self.name = name
        self.symbol = symbol
        self.system = system
        system.add_base(self)

    def __str__(self):
        return self.symbol

    def __hash__(self):
        return self._hash

    @property
    def factor(self):
        """
        The multiplication factor of a unit.
        For example, factor is 1000 for km. 
        """
        return 1

    @late
    def base(self):
        """
        The base represented as a list of tuples.
        Each tuple consists of an power and a unit.
        """
        return ((1,self),)

    def reference_string(self):
        return '{0}.base({1!r})'.format(self.system.reference_string(), self.quantity)


    def __reduce__(self):
        return (get_base_unit_with_name, (self.system, self.quantity,))


    def __eq__(self, other):
        if self is other:
            return True
        elif isinstance(other, base_unit):
            return NotImplemented
        else:
            return False


class no_system(object):
    ALL = {}

    @classmethod
    def set(cls, unit):
        cls.ALL[unit.name] = unit

    @classmethod
    def get(cls, name):
        return cls.ALL[name]

class none_unit(unit, metaclass=MultitonMetaClass):
    def __init__(self, name,  symbol):
        self.name = name
        self.symbol = symbol
        no_system.set(self)

    def __str__(self):
        return self.symbol


    def reference_string(self):
        return 'no_system.get({0!r})'.format(self.name)

    @late
    def factor(self):
        return 1

    @late
    def base(self):
        return ()

    def get_parts_with_power(self):
        return ()

    def is_none(self):
        return True


class zero_unit(none_unit):
    def __init__(self):
        none_unit.__init__(self,'zero', 'zero')

    def __str__(self):
        return self.symbol

    def is_zero(self):
        return True

    @late
    def base(self):
        return None

    def get_parts_with_power(self):
        return None


    def conversion_factor_from(self, x):
        if x.base is None:
            return 1.0
        else:
            return x.factor


class key_unit(none_unit):

    def iskey(self):
        return True

    @property
    def dtype(self):
        return numpy.dtype('uint64')




class nonnumeric_unit(unit):
    """
    nonnumeric_unit objects are  indivisable units 
    not connected to any system of units.
    
    nonnumeric_units cannot be used to
    derive new units from.
    
    nonnumeric_units have no physical meaning. 
    """
    def __init__(self, name, symbol):
        self.name = name
        self.symbol = symbol
        no_system.set(self)

    def __str__(self):
        return self.symbol

    def reference_string(self):
        return 'no_system.get({0!r})'.format(self.name)

    def __mul__(self, other):
        if other == 1:
            return self
        raise exceptions.AmuseException("Cannot derive other units from a non numeric unit")

    def __truediv__(self, other):
        raise exceptions.AmuseException("Cannot derive other units from a non numeric unit")

    def __rmul__(self, other):
        if other == 1:
            return self

        raise exceptions.AmuseException("Cannot derive other units from a non numeric unit")

    def __rtruediv__(self, other):
        if other == 1:
            return self

        raise exceptions.AmuseException("Cannot derive other units from a non numeric unit")

    def __pow__(self, other):
        if other == 1:
            return self

        raise exceptions.AmuseException("Cannot derive other units from a non numeric unit")

    def __div__(self, other):
        return self.__truediv__(other)

    def __rdiv__(self, other):
        return self.__rtruediv__(other)

    def is_non_numeric(self):
        return True

    @property
    def factor(self):
        return 1

    @property
    def base(self):
        return ((1,self),)

    def value_to_string(self, value):
        return None

    def is_valid_value(self, value):
        return False

    def __eq__(self, other):
        if self is other:
            return True
        elif isinstance(other, nonnumeric_unit):
            return NotImplemented
        else:
            return False

class string_unit(nonnumeric_unit):
    """
    String unit objects define quantities with a string value.
    These have no physical meaning, but are needed for some
    legacy codes. For example the path of a file.    
    """
    def __init__(self, name, symbol):
        nonnumeric_unit.__init__(self, name, symbol)

    def value_to_string(self, value):
        return '' if value is None else value

    def is_valid_value(self, value):
        return value is None or isinstance(value, str)

    @property
    def dtype(self):
        return numpy.dtype('S256')

class enumeration_unit(nonnumeric_unit):
    DEFINED={}
    """
    Enumeration unit objects define a fixed set of quantities.
    
    A quantity with a enumeration_unit can only have a
    value defined in the set of values of the enumeration_unit.
    
    :argument possible_values: A sequence or iterable with all 
        the possible values. If None the possible values are
        integers ranging from 0 to the length of the
        names_for_values argument
    :argument names_for_values: A sequence of strings defining a
        display name for each value. If None the names are the
        string vales of the values in the possible_values arguments
        
    Examples
    
    >>> my_unit = enumeration_unit('my_unit','my_unit', [1,2,5], ["star","gas","planet"])
    >>> 2 | my_unit
    quantity<2 - gas>
    >>> list(my_unit.quantities())
    [quantity<1 - star>, quantity<2 - gas>, quantity<5 - planet>]
    >>> 3 | my_unit
    Traceback (most recent call last):
        ...
    AmuseException: <3> is not a valid value for unit<my_unit>
    
    Or, with default values:
    
    >>> my_unit = enumeration_unit('my_unit','my_unit', None, ["star","gas","planet"])
    >>> 2 | my_unit
    quantity<2 - planet>
    >>> list(my_unit.quantities())
    [quantity<0 - star>, quantity<1 - gas>, quantity<2 - planet>]
    
    """
    def __init__(self, name, symbol, possible_values = None, names_for_values = None):
        nonnumeric_unit.__init__(self, name, symbol)

        self.possible_values = self._initial_list_of_possible_values(possible_values, names_for_values)
        self.names_for_values = self._initial_names_for_values(possible_values, names_for_values)
        if not len(self.possible_values) == len(self.names_for_values):
            raise exceptions.AmuseException("Must provide equal lenght list for values({0}) and names({1})".format(len(self.possible_values), len(self.names_for_values)))
        self.mapping_from_values_to_names = self._initial_mapping_from_values_to_names()
        self.DEFINED[name] = self

    def _initial_list_of_possible_values(self, possible_values, names_for_values):
        if possible_values is None:
            if names_for_values is None:
                raise exceptions.AmuseException("Must provide a list of values and / or a list of names for each value")
            return list(range(len(names_for_values)))
        else:
            return list(possible_values)


    def _initial_mapping_from_values_to_names(self):
        result = {}
        for value, name in zip(self.possible_values, self.names_for_values):
            result[value] = name
        return result


    def _initial_names_for_values(self, possible_values, names_for_values):
        if names_for_values is None:
            if possible_values is None:
                raise exceptions.AmuseException("Must provide a list of values and / or a list of names for each value")
            return [str(x) for x in possible_values]
        else:
            return list(names_for_values)

    def __hash__(self):
        return self._hash

    def is_valid_value(self, value):
        return value in self.mapping_from_values_to_names

    def value_to_string(self, value):
        return self.mapping_from_values_to_names[value]

    def quantities(self):
        for x in self.possible_values:
            yield x | self

    def __call__(self, string):
        index = self.names_for_values.index(string)
        if index > 0:
            return self.possible_values[index] | self
        else:
            raise exceptions.AmuseException("{0} is not a valid name for {1} enumeration type".format(string, self.name))

    @property
    def dtype(self):
        return numpy.dtype('int32')

    @classmethod
    def get(cls, name):
        try:
            return cls.DEFINED[name]
        except KeyError as ex:
            from amuse.units import nbody_system
            from amuse.units import si
            return cls.DEFINED[name]


    def __reduce__(self):
        return (get_enumeration_unit_with_name, (self.name,))





class named_unit(unit):
    """
    A named_unit object defines an alias for another
    unit. When printing a named_unit, the symbol
    is shown and not the unit parts. For all other
    operations the named_units works exactly like
    the aliased unit.
    
    :argument name: Long name or description of the unit
    :argument symbol: Short name to show when printing units
        or quantities
    :argument unit: The unit to alias
    
    >>> from amuse.units import si
    >>> 60.0 * si.s
    unit<60.0 * s>
    >>> minute = named_unit("minute","min", 60*si.s)
    >>> minute
    unit<min>
    >>> (20.0 | (60.0 * si.s)).as_quantity_in(minute)
    quantity<20.0 min>
    """
    def __init__(self, name, symbol, unit):
        self.name = name
        self.symbol = symbol
        self.local_unit = unit

    def __str__(self):
        return self.symbol


    def reference_string(self):
        return self.to_simple_form().reference_string()

    @late
    def factor(self):
        return self.local_unit.factor

    @late
    def base(self):
        return self.local_unit.base

    def is_none(self):
        return self.local_unit.is_none()

class derived_unit(unit, metaclass=MultitonMetaClass):
    """
    Abstract base class of derived units. New units
    can be derived from base_units. Each operation on
    a unit creates a new derived_unit.
    """
    pass


class factor_unit(derived_unit):
    """
    A factor_unit object defines a unit multiplied by
    a number. Do not call this method directly,
    factor_unit objects are supposed to be created by
    multiplying a number with a unit.
    
    :argument unit: The unit to derive from.
    :argument factor: The multiplication factor.
    
    >>> from amuse.units import si
    >>> minute = 60.0 * si.s
    >>> minute.as_quantity_in(si.s)
    quantity<60.0 s>
    >>> hour = 60.0 * minute
    >>> hour
    unit<60.0 * 60.0 * s>
    >>> hour.as_quantity_in(si.s)
    quantity<3600.0 s>
    
    """

    def __init__(self, factor, unit, name = None, symbol = None):
        self.name = name
        self.symbol = symbol
        self.local_factor = factor
        self.local_unit = unit

    def __str__(self):
        if self.symbol is None:
            return str(self.local_factor) + ' * ' + str(self.local_unit)
        return self.symbol + str(self.local_unit) 


    def reference_string(self):
        return '(' + str(self.local_factor) + ' * ' +  self.local_unit.reference_string() + ')'

    @late
    def factor(self):
        return self.local_factor * self.local_unit.factor

    @late
    def base(self):
        return self.local_unit.base


    def get_parts_with_power(self):
        local_unit_parts = self.local_unit.get_parts_with_power()
        result = []
        is_first = True
        for factor, power, unit in local_unit_parts:
            if is_first:
                factor *= self.local_factor
                is_first = False
            result.append( (factor, power, unit) )

        return result



class mul_unit(derived_unit):
    """
    A mul_unit object defines a unit multiplied by
    another unit. Do not call this method directly,
    mul_unit objects are supposed to be created by
    multiplying units.
    
    :argument left_hand: Left hand side of the multiplication.
    :argument right_hand: Right hand side of the multiplication.
    
    >>> from amuse.units import si
    >>> area = si.m * si.m
    >>> area
    unit<m * m>
    >>> hectare = (100 * si.m) * (100 * si.m)
    >>> hectare.as_quantity_in(area)
    quantity<10000.0 m * m>
    
    """
    def __init__(self, left_hand, right_hand):
        self.left_hand = left_hand
        self.right_hand = right_hand

    def __str__(self):
        return str(self.left_hand) + ' * ' + str(self.right_hand) 

    def reference_string(self):
        return '(' +  self.left_hand.reference_string() + ' * ' +  self.right_hand.reference_string() + ')'

    @late
    def factor(self):
        return self.left_hand.factor * self.right_hand.factor

    @late
    def base(self):
        return tuple(
            [
                x
                for x in [
                    (x[0] + x[1], x[2])
                    for x in self.combine_bases(self.left_hand.base, self.right_hand.base)
                ]
                if x[0] != 0
            ]
        )

    def get_parts_with_power(self):        
        lhs_parts = list(self.left_hand.get_parts_with_power())
        rhs_parts = list(self.right_hand.get_parts_with_power())

        result = []
        for lhs_factor, lhs_power, lhs_unit in lhs_parts: 
            rhs_index = 0
            found_match = False
            for rhs_factor, rhs_power, rhs_unit in rhs_parts:
                if lhs_unit is rhs_unit:
                    result.append( (lhs_factor * rhs_factor, lhs_power + rhs_power, lhs_unit,) )
                    found_match = True
                    del rhs_parts[rhs_index]
                    break
                rhs_index += 1
            if not found_match:
                result.append( (lhs_factor, lhs_power, lhs_unit,))

        for rhs_factor, rhs_power, rhs_unit in rhs_parts:
            result.append( (rhs_factor, rhs_power, rhs_unit,))
        return result




class pow_unit(derived_unit):
    """
    A pow_unit object defines a unit as
    another unit to a specified power. 
    
    Do not call this method directly,
    pow_unit objects are supposed to be created by
    taking powers of units.
    
    :argument power: Power of the unit
    :argument unit: The unit to derive from
    
    >>> from amuse.units import si
    >>> area = si.m**2
    >>> area
    unit<m**2>
    >>> area.as_quantity_in(si.m * si.m)
    quantity<1 m * m>
    >>> hectare = (100 * si.m) ** 2
    >>> hectare.as_quantity_in(area)
    quantity<10000.0 m**2>
    
    """
    def __init__(self, power, unit):
        self.power = power
        self.local_unit = unit

    def __str__(self):
        if isinstance(self.local_unit, derived_unit):
            return '(' + str(self.local_unit) + ')**' + str(self.power)
        else:
            return str(self.local_unit) + '**' + str(self.power)


    def reference_string(self):
        return '(' +  self.local_unit.reference_string() + '**' + str(self.power) + ')'

    @late
    def base(self):
        return tuple(
            [
                x
                for x in [
                    (x[0] * self.power, x[1])
                    for x in self.local_unit.base
                ]
                if x[0] != 0
            ]
        )

    @late
    def factor(self):
        return self.local_unit.factor ** self.power



    def get_parts_with_power(self):        

        result = []
        for factor, power, unit in self.local_unit.get_parts_with_power():
            result.append( (factor ** self.power, power * self.power, unit,))

        return result


class div_unit(derived_unit):
    """
    A div_unit object defines a unit multiplied by
    another unit. Do not call this method directly,
    div_unit objects are supposed to be created by
    dividing units.
    
    :argument left_hand: Left hand side of the multiplication.
    :argument right_hand: Right hand side of the multiplication.
    
    >>> from amuse.units import si
    >>> speed = si.m / si.s
    >>> speed
    unit<m / s>
    >>> speed_with_powers = si.m * si.s ** -1
    >>> speed.as_quantity_in(speed_with_powers)
    quantity<1 m * s**-1>
    
    """
    def __init__(self, left_hand, right_hand):
        self.left_hand = left_hand
        self.right_hand = right_hand

    def __str__(self):
        if isinstance(self.right_hand, derived_unit):
            return str(self.left_hand) + ' / (' + str(self.right_hand)+')'
        else:
            return str(self.left_hand) + ' / ' + str(self.right_hand)+''

    def reference_string(self):
        return '(' +  self.left_hand.reference_string() + '/' +  self.right_hand.reference_string() + ')'

    @late
    def factor(self):
        return  self.left_hand.factor * 1.0  / self.right_hand.factor

    @late
    def base(self):
        return tuple(
            [
                x
                for x in [
                    (x[0] - x[1], x[2])
                    for x in self.combine_bases(self.left_hand.base, self.right_hand.base)
                ]
                if x[0] != 0
            ]
        )


    def get_parts_with_power(self):        
        lhs_parts = list(self.left_hand.get_parts_with_power())
        rhs_parts = list(self.right_hand.get_parts_with_power())

        result = []
        for lhs_factor, lhs_power, lhs_unit in lhs_parts: 
            rhs_index = 0
            found_match = False
            for rhs_factor, rhs_power, rhs_unit in rhs_parts:
                if lhs_unit is rhs_unit:
                    result.append( (lhs_factor / rhs_factor, lhs_power - rhs_power, lhs_unit,) )
                    found_match = True
                    del rhs_parts[rhs_index]
                    break
                rhs_index += 1
            if not found_match:
                result.append( (lhs_factor, lhs_power, lhs_unit,))

        for rhs_factor, rhs_power, rhs_unit in rhs_parts:
            result.append( (1.0 / rhs_factor, -rhs_power, rhs_unit,))

        return result

class UnitException(exceptions.AmuseException):
    formatstring = "Unit exception: {0}"


class IncompatibleUnitsException(exceptions.AmuseException):
    formatstring = "Cannot express {1} in {0}, the units do not have the same bases"

    def __init__(self, *arguments):
        Exception.__init__(self)
        self.arguments = arguments

def get_system_with_name(name):
    return system.get(name)


def get_enumeration_unit_with_name(name):
    return enumeration_unit.get(name)

def get_base_unit_with_name(system, name):
    return system.base(name)



class UnitWithSpecificDtype(named_unit):

    def __init__(self, unit, dtype):
        self.specific_dtype = dtype
        symbol = str(unit) + "_" + str(dtype)
        named_unit.__init__(self, symbol, symbol, unit)

    @property
    def dtype(self):
        return self.specific_dtype

@memoize
def unit_with_specific_dtype(unit, dtype):
    if unit is None or dtype is None:
        return unit
    return UnitWithSpecificDtype(unit, dtype)

