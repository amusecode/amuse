"""
"""

from amuse.support.data import values

class attribute_def(object):
    def __init__(self, name, unit):
        self.name = name
        self.unit = unit
    
class Star(object):
    def __init__(self, id, **keyword_arguments):
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "attributes", set([]))
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        self.attributes.add(name_of_the_attribute)
        if isinstance(new_value_for_the_attribute, values.value):
            object.__setattr__(self, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise Exception("attribute "+name_of_the_attribute+" does not have a valid value, values must have a unit")
            
    def ensure_attributes_like(self, prototype_star):
        for x in prototype_star.attributes:
            if not x in self.attributes:
                setattr(self, x, getattr(prototype_star, x))
                
    def __str__(self):
        output = 'star<'
        output += str(self.id)
        output += '>'
        output += '\n'
        for x in self.attributes:
            output += x
            output += ': {'
            output += str(getattr(self,x))
            output += '}, '
            output += '\n'
        return output
    
    def to_si(self, convert_nbody):
        for x in self.attributes:
            value = getattr(self,x)
            setattr(self, x, convert_nbody.to_si(value))
