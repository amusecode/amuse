"""
"""
import numpy

class Quantity(object):
    
    def __init__(self, unit):
        self.unit = unit
        self.number = None
        
    def is_quantity(self):
        return True
        
    def is_scalar(self):
        return False
        
    def is_vector(self):
        return False
        
    def __repr__(self):
        return 'quantity<'+str(self)+'>'

    def __div__(self, other):
        return new_quantity(self.number / other.number , (self.unit / other.unit).to_simple_form())
        
    def __add__(self, other):
        other_in_my_units = other.in_(self.unit)
        return new_quantity(self.number + other_in_my_units.number , self.unit)
        
    def __sub__(self, other):
        other_in_my_units = other.in_(self.unit)
        return new_quantity(self.number - other_in_my_units.number , self.unit)

    def __mul__(self, other):
        return  new_quantity(self.number * other.number , (self.unit * other.unit).to_simple_form())
        
            
    def in_(self, another_unit):
        value_of_unit_in_another_unit = self.unit.in_(another_unit)
        return new_quantity(self.number * value_of_unit_in_another_unit.number, another_unit)

    def value_in(self, another_unit):
        value_of_unit_in_another_unit = self.unit.in_(another_unit)
        return self.number * value_of_unit_in_another_unit.number

class ScalarQuantity(Quantity):
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
        other_in_my_units = other.in_(self.unit)
        return self.number < other_in_my_units.number
        
    def __gt__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number > other_in_my_units.number
        
    def __eq__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number == other_in_my_units.number
        
    def __neq__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number != other_in_my_units.number

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)
            
            
class VectorQuantity(Quantity):
    def __init__(self, array, unit):
        Quantity.__init__(self, unit)
        self.x = array[0]
        self.y = array[1]
        self.z = array[2]
          
    def is_vector(self):
        return True
        
    @property
    def number(self):
        return numpy.array([self.x, self.y, self.z])
        
    def __str__(self):
        unit_str = str(self.unit)
        if unit_str:
            return '[' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ']' + ' | ' + unit_str
        else:
            return '[' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ']'
            

    def __lt__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number < other_in_my_units.number
        
    def __gt__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number > other_in_my_units.number
        
    def __eq__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number == other_in_my_units.number
        
    def __neq__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number != other_in_my_units.number        
                 
        
def new_quantity(value, unit):
    if isinstance(value, list):
       return VectorQuantity(value, unit)
    if isinstance(value, numpy.ndarray):
       return VectorQuantity(value, unit)
    return ScalarQuantity(value, unit)
       
