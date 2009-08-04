"""
"""

class value(object):
    def __init__(self, number, unit):
        self.number = number
        self.unit = unit
        
    def __str__(self):
        unit_str = str(self.unit)
        if unit_str:
            return str(self.number) + ' ' + unit_str
        else:
            return str(self.number)
            
    def __repr__(self):
        return 'value<'+str(self)+'>'
        
    def __div__(self, other):
        return (self.unit / other.unit)(self.number / other.number)
    
    def __add__(self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number + other_in_my_units.number | self.unit
        
    def __lt__ (self, other):
        other_in_my_units = other.in_(self.unit)
        return self.number < other_in_my_units.number
            
    def in_(self, another_unit):
        value_of_unit_in_another_unit = self.unit.in_(another_unit)
        return another_unit(self.number * value_of_unit_in_another_unit.number)
        
class values(object):
    def __init__(self, numbers, unit):
        self.numbers = numbers
        self.unit = unit
    def __getitem__(self, index):
        return self.unit(self.numbers[index]) 
       
