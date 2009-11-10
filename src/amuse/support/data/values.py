"""
"""
import numpy

class Quantity(object):
    
    def __init__(self, unit):
        self.unit = unit
        
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
        self._number = numpy.array(array)
        
    def is_vector(self):
        return True
    
    
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
        quantity_in_my_units = quantity.in_(self.unit)
        self._number[index] = quantity_in_my_units.number
        
    @property
    def number(self):
        return self._number
        
    @property
    def x(self):
        """The x axis component of a 3 dimensional vector.
        This is equavalent to the first item in an array.
        
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
        This is equavalent to the first item in an array.
        
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
        This is equavalent to the first item in an array.
        
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
            
    def indices(self):
        for x in len(self._number):
            yield x
                 
        
def new_quantity(value, unit):
    if isinstance(value, list):
       return VectorQuantity(value, unit)
    if isinstance(value, numpy.ndarray):
       return VectorQuantity(value, unit)
    return ScalarQuantity(value, unit)
       
