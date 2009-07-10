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
    def __str__(self):
        output = 'star<'
        output += str(self.id)
        output += '>'
        for x in self.attributes:
            output += x
            output += ': '
            output += str(getattr(self,x))
            output += ', '
        return output