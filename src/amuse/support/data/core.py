"""
"""

from amuse.support.data import values
from amuse.support.units import si

import numpy

class TemporalAttribute(object):
    def __init__(self, name):
        self.values = []
        self.name = name
    
    def get_times(self):
        for time, value in self.values:
            yield time
            
    def set_value_at_time(self, time, value):
        self.values.append((time, value))
    
    def get_value_at_time(self, requested_time = None):
        if requested_time is None:
            return self.value()
            
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
        return self.value().value_in(units)
        
    def __str__(self):
        return str(self.time()) + " - " + str(self.value())
    
class AttributeList(object):
    
    def __init__(self, particles, attributes, lists_of_values, units):
        d = {}
        for index, particle in enumerate(particles):
            d[id(particle)] = index
        self.mapping_from_particle_to_index = d
        
        self.mapping_from_attribute_to_values_and_unit = {}
        for attribute, values, unit in zip(attributes, lists_of_values, units):
            self.mapping_from_attribute_to_values_and_unit[attribute] = (unit,values)
        
        self.particles = particles
        
    def get_value_of(self, particle, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            return None
        unit,values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        index = self.mapping_from_particle_to_index.get(id(particle))
        if index is None:
            return None
        else:
            return values[index] | unit
            
    def iter_values_of(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            return
            
        unit,values = self.mapping_from_attribute_to_values_and_unit[attribute]
        particles = self.particles
        
        for index in range(len(self.particles)):
            yield particles[i], (values[i] | unit)
    
    
    def get_indices_of(self, particles):
        mapping_from_particle_to_index = self.mapping_from_particle_to_index 
        result = numpy.zeros(len(particles),dtype='int32')
        #result = [mapping_from_particle_to_index[id(particle)] for particle in particles]
        
        index = 0
        for index, particle in enumerate(particles):
            result[index] = mapping_from_particle_to_index[id(particle)]
            index += 1
        return result
        
    
    
    def get_values_of_particles_in_units(self, particles, attributes, target_units):
        indices = self.get_indices_of(particles)
        results = []
        for attribute, target_unit in zip(attributes, target_units):
             unit, values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_unit_in_target_unit = unit.value_in(target_unit )
             selected_values = values.take(indices)
             if value_of_unit_in_target_unit != 1.0:
                selected_values *= value_of_unit_in_target_unit
             results.append(selected_values)
        
        return results
    
    
    def set_values_of_particles_in_units(self, particles, attributes, list_of_values_to_set, source_units):
        indices = self.get_indices_of(particles)
        
        results = []
        for attribute, values_to_set, source_unit in zip(attributes, list_of_values_to_set, source_units):
             self.set_values_of_indexed_attribute_in_units(indices, attribute, values_to_set, source_unit)
                 
        return results
        
    def set_values_of_indexed_attribute_in_units(self, indices, attribute, values_to_set, source_unit):
        if attribute in self.mapping_from_attribute_to_values_and_unit:
             unit, values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_source_unit_in_list_unit = source_unit.value_in(unit)
             selected_values = values_to_set
             if value_of_source_unit_in_list_unit != 1.0:
                 selected_values *= value_of_source_unit_in_list_unit
             values.put(indices, selected_values)
        else:
             values = numpy.empty(len(self.particles))
             unit = source_unit
             self.mapping_from_attribute_to_values_and_unit[attribute] = (unit, values)
             selected_values = values_to_set
             values.put(indices, selected_values)
             
    def merge_into(self, others):
        source_attributes = []
        source_units = []
        source_valeus = []
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            source_unit, source_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            source_attributes.append(attribute)
            source_values.append(source_values)
            source_units.append(source_units)
            
                
        other.set_values_of_particles_in_units(self.particles, source_attributes, source_values, source_units)
        
    def remove_particles(self, particles):
        indices = self.get_indices_of(particles)
        
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            unit, values = self.mapping_from_attribute_to_values_and_unit[attribute]
            values.delete(indices)
        
        
        

        
        
    
    
    
class Particles(object):
    """A set of particle objects"""
    def __init__(self, size = 0):
        self.particles = [Particle(i+1, self) for i in range(size)]
        #selt.attributelist = AttributeList(
    def __iter__(self):
        return iter(self.particles)

    def __getitem__(self, index):
        return self.particles[index]
        
    @property
    def mass(self):
        result = []
        for x in self:
            result.append((x.id, x.mass))
        return result
    
    @property
    def ids(self):
        for x in self.particles:
            yield x.id

    def get_values_of_attribute(self, attribute, time = None):
        result = []
        for x in self:
            result.append(getattr(x, attribute).get_value_at_time(time))
        return result
        
    def set_values_of_attribute(self, attribute, time, values):
        for particle, value in map(None, self, values):
            getattr(particle, attribute).set_value_at_time(time, value)
        
    def ids_for_module_with_id(self, module_id):
        for x in self.particles:
            if module_id in x._module_ids_to_index:
                yield  x._module_ids_to_index[module_id][1]
            else:
                pass
    
    
    def get_values_of_attribute_for_module_with_id(self, module_id, attribute, time = None):
        result = []
        for x in self:
            if module_id in x._module_ids_to_index:
                result.append(getattr(x, attribute).get_value_at_time(time))
            else:
                pass
        return result
        
    
    
    def set_values_of_attribute_for_module_with_id(self, module_id, attribute, values, time = None, times = None):
        result = []
        index = 0
        for x in self:
            if module_id in x._module_ids_to_index:
                if times is not None:
                    time = times[index]
                if not hasattr(x, attribute):
                    setattr(x, attribute, values[index])
                getattr(x, attribute).set_value_at_time(time, values[index])
                index += 1
            else:
                pass
        return result
        
class Stars(Particles):
    
    def center_of_mass(self):
        sum_mass = 0.0 | si.kg
        sum_massposition = [0.0, 0.0, 0.0] | si.kg * si.m
        for star in self:
            sum_mass += star.mass.value()
            sum_massposition += star.position.value() * star.mass.value()
        return sum_massposition / sum_mass
        
class Measurement(object):
    def __init__(self, timestamp,  attributes, units,  ids, values=None):
        self.timestamp = timestamp
        self.ids = ids
        self.attributes = attributes
        if values == None:
            self.values = zeros((len(attributes), len(ids)))
        else:
            self.values = values
        self.units = units
        
        if self.values.ndim != 2:
            raise Exception("values must be a 2 dimensional array")
        if self.values.shape[0] != len(self.attributes):
            raise Exception("there must be the same number of columns in the values array as there are attributes")
        if len(self.attributes) != len(self.units):
            raise Exception("there must be the same number of attributes as units")
        if self.values.shape[1] != len(ids):
            raise Exception("there must be the same number of rows in the values array as there are ids")
            
        self.ids_to_rownumber = {}
        for rownumber, id in enumerate(self.ids):
            self.ids_to_rownumber[id] = rownumber
            
        self.attribute_to_colnumber = {}
        for colnumber, attribute in enumerate(self.attributes):
            self.attribute_to_colnumber[attribute] = colnumber
            
    def __str__(self):
        rows = []
        columns = map(lambda x : str(x), self.attributes)
        columns.insert(0,'id')
        rows.append(columns)
        
        columns = map(lambda x : str(x), self.units)
        columns.insert(0,'-')
        rows.append(columns)
        rows.append(map(lambda x : '========', range(len(self.units)+1)))
        for i in range(self.values.shape[1]):
            row = [str(self.ids[i])]
            for j in range(self.values.shape[0]):
                row.append(str(self.values[j,i]))
            rows.append(row)
        lines = map(lambda  x : '\t'.join(x), rows)
        return '\n'.join(lines)
        
    def get_value(self, attribute, id):
        rownumber = self.ids_to_rownumber[id]
        colnumber = self.attribute_to_colnumber[attribute]
        return self.values[colnumber][rownumber] | self.units[colnumber]
            
            
class Particle(object):
    """A physical object or a physical region simulated as a 
    physical object (cload particle).
    
    All attributes defined on a particle are specific for 
    that particle (for example mass or position). A particle contains 
    a set of attributes, some attributes are *generic* and applicaple
    for multiple modules. Other attributes are *specific* and are 
    only applicable for a single module.
    """
    
    def __init__(self, id = -1, set = None, **keyword_arguments):
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "attributes", {})
        object.__setattr__(self, "_module_ids_to_index", {})
        if not set is None:
            object.__setattr__(self, "set", set)
            
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if not name_of_the_attribute in self.attributes:
            self.attributes[name_of_the_attribute] = TemporalAttribute(name_of_the_attribute)
        if isinstance(new_value_for_the_attribute, values.Quantity):
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
            
    def set_value_of_attribute(self, attribute, value, time = None):
        getattr(self, attribute).set_value_at_time(time, values[index])
            
class Star(Particle):
    pass
