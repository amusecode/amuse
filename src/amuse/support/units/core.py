from amuse.support.data import values

import numpy

class system(object):
    def __init__(self, name):
        self.name = name
        self.bases = []
    
    def add_base(self, unit):
        unit.system = self
        unit.index = len(self.bases)
        self.bases.append(unit)
        
class unit(object):
    """
    Abstract base class for unit objects.
        
    Two classes of units are defined:

    base units
        The base units in a given system of units. For SI, these
        are meter, kilogram, second, ampere, kelvin, mole and 
        candele. See the si module :mod:`amuse.support.units.si`
    derived units
        Derived units are created by dividing or multiplying
        with a number or with another unit. For example, 
        to get a velocity unit we can devine vel = 1000 * m / s

    Units can also be named, by creating a named unit.
    """
    
    def __mul__(self, other):
       
        if isinstance(other, unit):
            return mul_unit(self, other)
        else:
            return factor_unit(other, self)
        
    def __div__(self, other):
        if isinstance(other, unit):
            return div_unit(self, other)
        else:
            return factor_unit(1.0 / other, self)

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
        
        >>> from amuse.support.units import units
        >>> 100 | units.kg
        quantity<100 kg>
        """
        return self.new_quantity(value)
        
    def __rdiv__(self, other):
        return factor_unit(other, pow_unit(-1,self))
        
    def __pow__(self, other):
        if other == 1:
            return self
        else:
            return pow_unit(other, self)
        
    def __call__(self, x):
        return self.new_quantity(x)
        
    
    def new_quantity(self, value):
        """Create a new Quantity object.
    
        :argument value: numeric value of the quantity, can be 
            a number or a sequence (list or ndarray)
        :returns: new ScalarQuantity or VectorQuantity object 
            with this unit
        """
        return values.new_quantity(value, self)
        
    def to_simple_form(self):
        """Convert unit to a form with only one factor and powers
        
        :result: Unit with only a factor and power terms
        
        >>> from amuse.support.units import units
        >>> N = (units.m * units.kg) / (units.s * units.s)
        >>> N
        unit<m * kg / s * s>
        >>> J = N * units.m
        >>> J
        unit<m * kg / s * s * m>
        >>> J.to_simple_form()
        unit<m**2 * kg * s**-2>
        """
        
        if not self.base:
            return none_unit('no unit', '')
            
        result = self.factor
        for n, base in self.base:
            result =  result * (base ** n)
        
        return result
    
    def are_bases_equal(self, other):
        for n1, unit1 in sort(self.base, key=lambda x: x[1].index):
            found = False
            for n2, unit2 in other.base:
                if unit1 == unit2:
                    if not n2 == n1:
                        return False
                    found = True
                    break;
            if not found:
                return False
        return True
                        
    def conversion_factor_from(self, x):
        if self.base == x.base:
            this_factor = self.factor * 1.0
            other_factor = x.factor
            return this_factor / other_factor
        else:
            raise Exception("Cannot expres: " + str(x) + " in " + str(self))
      
    def in_(self, x):
        return self.as_quantity_in(x)
    
    def as_quantity_in(self, unit):
        """Express this unit as a quantity in the given unit
        
        :argument unit: The unit to express this unit in
        :result: A Quantity object
        
        Examples
        
        >>> from amuse.support.units import units
        >>> ton = 1000 * units.kg
        >>> ton.as_quantity_in(units.kg)
        quantity<1000.0 kg>
        """
        if isinstance(unit, values.Quantity):
            raise Exception("Cannot expres a unit in a quantity")
        else:
            factor = self.conversion_factor_from(unit)
            return values.new_quantity(factor, unit)
            
    def value_in(self, unit):
        """
        Return a numeric value of this unit in the given unit.
        Works only when the units are compatible, i.e. from
        tonnage to kg's.
        
        A number is returned without any unit information.
        
        :argument unit: wanted unit of the value
        :returns: number in the given unit
        
        >>> from amuse.support.units import units
        >>> x = units.km
        >>> x.value_in(units.m)
        1000.0
        
        """
        return self.conversion_factor_from(unit)
        
    def __repr__(self):
        return 'unit<'+str(self)+'>'
        
    def combine_bases(self, base1, base2):
        result = [None] * 7
        for n1, unit1 in base1:
            found = False
            for n2, unit2 in base2:
                if unit1 == unit2:
                    base2 = filter(lambda x : x[1] != unit1, base2)
                    found = True
                    result[unit1.index] = (n1, n2, unit1)
                    break
            if not found:
                result[unit1.index] =  n1, 0, unit1
        for n2, unit2 in base2:
                result[unit2.index] = 0, n2, unit2
        for x in result:
            if not x is None:
                yield x
                
    def has_same_base_as(self, other):
        """Detrmine if the base of other is the same as the
        base of self.
        
        :argument other: unit to compare base to
        :result: True, if bases are compatiple.
        
        >>> from amuse.support.units import units
        >>> mps = units.m / units.s
        >>> kph = units.km / units.hour
        >>> mps.has_same_base_as(kph)
        True
        >>> mps.has_same_base_as(units.km)
        False
        
        """
        return other.base == self.base
            
        
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
    def __init__(self, quantiy, name, symbol, system):
        self.quantity = quantiy
        self.name = name
        self.symbol = symbol
        system.add_base(self)
        
    def __str__(self):
        return self.symbol
    
    @property
    def factor(self):
        """
        The multiplication factor of a unit.
        For example, factor is 1000 for km. 
        """
        return 1
        
    @property
    def base(self):
        """
        The base represented as a list of tuples.
        Each tuple consists of an power and a unit.
        """
        return ((1,self),)
        
class none_unit(unit):
    def __init__(self, name, symbol):
        self.name = name
        self.symbol = symbol
        
    def __str__(self):
        return self.symbol
    
    @property
    def factor(self):
        return 1.0
        
    @property
    def base(self):
        return ()
        
    
class named_unit(unit):
    def __init__(self, name, symbol, unit):
        self.name = name
        self.symbol = symbol
        self.unit = unit
        
    def __str__(self):
        return self.symbol
    
    @property
    def factor(self):
        return self.unit.factor
        
    @property
    def base(self):
        return self.unit.base
        

class derived_unit(unit):
    pass
    
class factor_unit(derived_unit):
    def __init__(self, factor , unit, name = None, symbol = None):
        self.name = name
        self.symbol = symbol
        self.local_factor = factor
        self.unit = unit
        
    def __str__(self):
        if self.symbol is None:
            return str(self.local_factor) + ' * ' + str(self.unit)
        return self.symbol + str(self.unit) 
    
    @property
    def factor(self):
        return self.local_factor * self.unit.factor
        
    @property
    def base(self):
        return self.unit.base
        
class mul_unit(derived_unit):
    
    def __init__(self, left_hand , right_hand):
        self.left_hand = left_hand
        self.right_hand = right_hand
        
    def __str__(self):
        return str(self.left_hand) + ' * ' + str(self.right_hand) 
        
    @property
    def factor(self):
        return self.left_hand.factor * self.right_hand.factor
   
    @property
    def base(self):
        return tuple(
            filter(lambda x: x[0] != 0,
            map(lambda x: (x[0] + x[1], x[2]),self.combine_bases(self.left_hand.base, self.right_hand.base))))
        
class pow_unit(derived_unit):
    
    def __init__(self, power , unit):
        self.power = power
        self.unit = unit
        
    def __str__(self):
        return str(self.unit) + '**' + str(self.power)
        
    @property
    def base(self):
        return tuple(
            filter(lambda x: x[0] != 0,
            map(lambda x : (x[0] * self.power, x[1]), self.unit.base)))
        
    @property
    def factor(self):
        return self.unit.factor ** self.power
        
class div_unit(derived_unit):
    
    def __init__(self, left_hand , right_hand):
        self.left_hand = left_hand
        self.right_hand = right_hand
        
    def __str__(self):
        return str(self.left_hand) + ' / ' + str(self.right_hand)+''
        
    @property
    def factor(self):
        return  self.left_hand.factor * 1.0  / self.right_hand.factor
        
    @property
    def base(self):
        return tuple(
                    filter(lambda x: x[0] != 0,
                    map(lambda x: (x[0] - x[1], x[2]),
                        self.combine_bases(self.left_hand.base, self.right_hand.base))))
        
def k(unit):
    return factor_unit(1000, unit, 'kilo','k')           
    

    
