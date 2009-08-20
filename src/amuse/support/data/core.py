"""
"""

from amuse.support.data import values
from amuse.support.units import si

class TemporalAttribute(object):
    def __init__(self, name):
        self.values = []
        self.name = name
    
    def get_times(self):
        for time, value in self.values:
            yield time
            
    def set_value_at_time(self, time, value):
        self.values.append((time, value))
    
    def get_value_at_time(self, requested_time):
        min = 0
        max = len(self.values) - 1
        while True:
            index_of_midpoint = (max + min) / 2 
            time, value = self.values[index_of_midpoint]
            requested_time_is_after_time_of_index = requested_time > time
            if max - min == 1:
                if requested_time_is_after_time_of_index:
                    return self.values[max]
                else:
                    return self.values[min]
                    
            if requested_time_is_after_time_of_index:
                min = index_of_midpoint
            else:
                max = index_of_midpoint
            
    
    def value(self):
        return self.values[-1][1]
    
    def time(self):
        return self.values[-1][0]
        
    def to_number_in(self, units):
        return self.value().in_(units).number
        
    def __str__(self):
        return str(self.time()) + " - " + str(self.value())
    
class Particles(object):
    """A set of particle objects"""
    def __init__(self, particles):
        self.particles = particles
        
    def __iter__(self):
        return iter(self.particles)

    @property
    def mass(self):
        result = []
        for x in self:
            result.append((x.id, x.mass))
        return result
   
class Particle(object):
    """A physical object or a physical region simulated as a 
    physical object (cload particle).
    
    All attributes defined on a particle are specific for 
    that particle (for example mass or position). A particle contains 
    a set of attributes, some attributes are *generic* and applicaple
    for multiple modules. Other attributes are *specific* and are 
    only applicable for a single module.
    """
    
    def __init__(self, id = -1, **keyword_arguments):
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "attributes", {})
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if not name_of_the_attribute in self.attributes:
            self.attributes[name_of_the_attribute] = TemporalAttribute(name_of_the_attribute)
        if isinstance(new_value_for_the_attribute, values.value):
            self.attributes[name_of_the_attribute].set_value_at_time( 0 | si.s ,new_value_for_the_attribute)
        else:
            raise Exception("attribute "+name_of_the_attribute+" does not have a valid value, values must have a unit")
    
    def __getattr__(self, name):
         return self.attributes[name] 
         
    def ensure_attributes_like(self, prototype_star):
        for name, attribute in prototype_star.attributes.iteritems():
            if not name in self.attributes:
                setattr(self, name, attribute.value())
                
    def __str__(self):
        output = 'Particle '
        output += str(self.id)
        output += ''
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
