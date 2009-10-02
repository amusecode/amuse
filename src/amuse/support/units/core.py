
"""
Unit conversion support
=======================

Definition of the 7 base units

"""
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
    
    def __ror__(self, x):
        return values.new_quantity(x, self)
        
    def __rdiv__(self, other):
        return factor_unit(other, pow_unit(-1,self))
        
    def __pow__(self, other):
        if other == 1:
            return self
        else:
            return pow_unit(other, self)
        
    def __call__(self, x):
        return values.new_quantity(x, self)
        
    def to_simple_form(self):
        result = self.factor
        for n, base in self.base:
            result =  result * (base ** n)
        return result
    
    def are_bases_equal(self, other):
        for n1, unit1 in self.base:
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
        if self.are_bases_equal(x):
            this_factor = self.factor * 1.0
            other_factor = x.factor
            return this_factor / other_factor
        else:
            print x.base, self.base
            raise Exception("Cannot expres: " + str(x) + " in " + str(self))
            
    def in_(self, x):
        if isinstance(x, values.Quantity):
            print "bla"
        else:
            factor = self.conversion_factor_from(x)
            return values.new_quantity(factor, x)
            
    def __repr__(self):
        return 'unit<'+str(self)+'>'
        
    def combine_bases(self, base1, base2):
        result = []
        for n1, unit1 in base1:
            found = False
            for n2, unit2 in base2:
                if unit1 == unit2:
                    base2 = filter(lambda x : x[1] != unit1, base2)
                    found = True
                    yield n1, n2, unit1
                    break
            if not found:
                yield n1, 0, unit1
        for n2, unit2 in base2:
                yield 0, n2, unit2
                
    def has_same_base_as(self, other):
        return other.base == self.base
            
        
class base_unit(unit):
    def __init__(self, quantiy, name, symbol, system):
        self.quantity = quantiy
        self.name = name
        self.symbol = symbol
        system.add_base(self)
    def __str__(self):
        return self.symbol
    @property
    def factor(self):
        return 1
    @property
    def base(self):
        return ((1,self),)
class none_unit(unit):
    def __init__(self, name, symbol):
        self.name = name
        self.symbol = symbol
    def __str__(self):
        return self.symbol
    @property
    def factor(self):
        return 1
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
        
class factor_unit(unit):
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
        
class mul_unit(unit):
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
        return tuple(map(lambda x: (x[0] + x[1], x[2]),self.combine_bases(self.left_hand.base, self.right_hand.base)))
        
class pow_unit(unit):
    def __init__(self, power , unit):
        self.power = power
        self.unit = unit
    def __str__(self):
        return str(self.unit) + '**' + str(self.power)
        
    @property
    def base(self):
        return tuple(map(lambda x : (x[0] * self.power, x[1]), self.unit.base))
        
    @property
    def factor(self):
        return self.unit.factor ** self.power
        
class div_unit(unit):
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
    

    
