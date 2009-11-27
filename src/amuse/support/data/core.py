"""
"""

from amuse.support.data import values
from amuse.support.units import si

import numpy

class AttributeValues(object):
    __slots__ = ["attribute", "values", "unit", "model_times"]
    
    def __init__(self, attribute, unit, values = None,  model_times = None, length = None):
        self.attribute = attribute
        self.unit = unit
        self.model_times = model_times
        if values is None:
            self.values = numpy.zeros(length)
        else:
            self.values = values
        
    def copy(self):
        return AttributeValues(self.attribute, self.unit, self.values.copy(), self.model_times)
        
class AttributeList(object):
    
    def __init__(self, particles, attributes = [], lists_of_values = [], units = [], model_times = None):
        
        if len(lists_of_values) != len(attributes):
            raise Exception(
                "you need to provide the same number of value list as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(lists_of_values)
                )
            )
        if len(attributes) != len(units):
             raise Exception(
                "you need to provide the same number of value list as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(units)
                )
            )
        if len(lists_of_values) > 0 and len(particles) != len(lists_of_values[0]):
            raise Exception(
                "you need to provide the same number of values as particles, found {0} values and {1} particles".format(
                    len(lists_of_values[0]), len(particles)
                )
            )
        
        
        model_times = self._convert_model_times(model_times, particles)
        
        self.mapping_from_attribute_to_values_and_unit = {}
        for attribute, values, unit in zip(attributes, lists_of_values, units):
            self.mapping_from_attribute_to_values_and_unit[attribute] = AttributeValues(
                attribute,
                unit,
                values,
                model_times
            )
        
        self.particles = numpy.array(particles)
        self.reindex()
        
    def copy(self):
        copy = AttributeList([])
        copy.mapping_from_particle_to_index = self.mapping_from_particle_to_index.copy()
        copy.particles = self.particles.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_values_and_unit.iteritems():
            copy.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            raise AttributeError("particle does not have a "+attribute)
        
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        index = self.mapping_from_particle_to_index[id(particle)]
        
        return attribute_values.values[index] | attribute_values.unit
        
    def get_values_of_attributes_of_particle(self, particle, attributes):
        
        index = self.mapping_from_particle_to_index[id(particle)]
        result = []
        for attribute in attribute:
            if not attribute in self.mapping_from_attribute_to_values_and_unit:
                raise AttributeError("particle does not have a "+attribute)
            
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            
            result.append(attribute_values.values[index] | attribute_values.unit)
        return result
        
    
    def set_value_of(self, particle, attribute, quantity):
        index = self.mapping_from_particle_to_index[id(particle)]
            
        if index is None:
            raise Exception("unknown particle " + particle)
            
        attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
             
        value_to_set = quantity.value_in(attribute_values.unit)
        attribute_values.values[index] = value_to_set
        
            
    def iter_values_of_particle(self, particle):
        index = self.mapping_from_particle_to_index[id(particle)]
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            yield attribute, (attribute_values.values[index] | attribute_values.unit)
    
    
            
    def iter_values_of(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            return
            
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        values = attribute_values.values
        unit = attribute_values.unit
        particles = self.particles
        
        for index in range(len(self.particles)):
            yield particles[i], (values[i] | unit)
            
    def iter_particles_as_views(self):
        class ParticleView(object):
            __slots__=['index']
            pass
        
        class AttributeViewProperty(object):
            __slots__ = ['attribute', 'values','unit']
            def __init__(self, list, attribute):
                self.attribute = attribute
                self.values = list.mapping_from_attribute_to_values_and_unit[attribute].values
                self.unit = list.mapping_from_attribute_to_values_and_unit[attribute].unit
            
            def __get__(self, instance, owner):
                return  values.ScalarQuantity(self.values[instance.index] , self.unit)
                
            def __set__(self, instance, value):
                self.values[instance.index] = value.value_in(self.unit)
                
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            setattr(
                ParticleView, 
                attribute, 
                AttributeViewProperty(self, attribute)
            )
            
        index = 0
        for index in range(len(self.particles)):
            p = ParticleView()
            p.index = index
            yield p
            index += 1
            
            
    def iter_particles(self):
        for particle in self.particles:
            yield particle
    
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
            
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_unit_in_target_unit = attribute_values.unit.value_in(target_unit )
             selected_values = attribute_values.values.take(indices)
             if value_of_unit_in_target_unit != 1.0:
                selected_values *= value_of_unit_in_target_unit
             results.append(selected_values)
        
        return results
        
    def get_values_of_all_particles_in_units(self, attributes, target_units):
        results = []
        for attribute, target_unit in zip(attributes, target_units):
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_unit_in_target_unit = attribute_values.unit.value_in(target_unit )
             selected_values = attribute_values.values
             if value_of_unit_in_target_unit != 1.0:
                selected_values = selected_values * value_of_unit_in_target_unit
             results.append(selected_values)
        
        return results
    
    def _convert_model_times(self, model_times, particles):
        if not model_times is None and isinstance(model_times, values.ScalarQuantity):
            return numpy.linspace(model_times.number, model_times.number, len(self.particles)) | model_times.unit
        else:
            return model_times
    
    def set_values_of_particles_in_units(self, particles, attributes, list_of_values_to_set, source_units, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, particles)
        
        previous_model_times = None
        results = []
        for attribute, values_to_set, source_unit in zip(attributes, list_of_values_to_set, source_units):
            selected_values = values_to_set
            
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                   attribute,
                   source_unit,
                   length = len(self.particles)
                )
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                 
            value_of_source_unit_in_list_unit = source_unit.value_in(attribute_values.unit)
            if value_of_source_unit_in_list_unit != 1.0:
                selected_values *= value_of_source_unit_in_list_unit 
             
            attribute_values.values.put(indices, selected_values)
            if not model_times is None:
                if not previous_model_times is attribute_values.model_times:
                    attribute_values.model_times.put(indices, model_times)
                    previous_model_times = attribute_values.model_times
            
        return results
        
             
    def merge_into(self, others):
        source_attributes = []
        source_units = []
        source_valeus = []
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            source_attributes.append(attribute_values.attribute)
            source_values.append(attribute_values.values)
            source_units.append(attribute_values.unit)
            
                
        other.set_values_of_particles_in_units(self.particles, source_attributes, source_values, source_units)
        
    def remove_particles(self, particles):
        indices = self.get_indices_of(particles)
        
        mapping_from_attribute_to_values_and_unit = self.mapping_from_attribute_to_values_and_unit.copy()
        for attribute, attribute_values in mapping_from_attribute_to_values_and_unit.iteritems():
            attribute_values.values = numpy.delete(attribute_values.values,indices)
        
        self.particles = numpy.delete(self.particles,indices)
        self.reindex()
        
    def reindex(self):
        d = {}
        index = 0
        for particle in self.particles:
            d[id(particle)] = index
            index += 1
          
        self.mapping_from_particle_to_index = d
        
    def get_attribute_as_vector(self, particle, names_of_the_scalar_attributes):
        index = self.mapping_from_particle_to_index[id(particle)]
        vector = []
        for i, attribute in enumerate(names_of_the_scalar_attributes):
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            vector.append(attribute_values.values[index])
            unit_of_the_vector = attribute_values.unit
        return vector | unit_of_the_vector
        
    def get_attributevalues_for(self, attribute, unit):
        if attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        else:
            attribute_values = AttributeValues(
                attribute,
                unit,
                length = len(self.particles)
            )
            self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
        return attribute_values
        
    def set_attribute_as_vector(self, particle, names_of_the_scalar_attributes, quantity):
        index = self.mapping_from_particle_to_index[id(particle)]
        vector = []
        for i, attribute  in enumerate(names_of_the_scalar_attributes):
            attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
        
            value_to_set = quantity[i].value_in(attribute_values.unit)
            attribute_values.values[index] = value_to_set
            
    def attributes(self):
        return set(self.mapping_from_attribute_to_values_and_unit.keys())
    
    def __str__(self):
        attributes = sorted(self.attributes)
        
        columns = map(lambda x : [x], attributes)
        columns.insert(0,['id'])
        
        for index, attribute in enumerate(attributes):
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            column = columns[index + 1]
            column.append(str(attribute_values.unit))
            column.append('========')
            if len(attribute_values.values) > 40:
                values_to_show = list(attribute_values.values[:20])
                values_to_show.append(attribute_values.values[-20:])
            else:
                values_to_show = attribute_values.values
            
            for value in values_to_show:
                column.append(str(value))
            column.append('========')
            
        column = columns[0]
        column.append("-")
        column.append('========')
        if len(self.particles) > 40:
            values_to_show = list(self.particles[:20])
            values_to_show.append(self.particles[-20:])
        else:
            values_to_show = self.particles
                    
        for value in values_to_show:
            column.append(str(id(value)))
            
        column.append('========')
            
        rows = []
        for i in range(len(columns[0])):
            row = [x[i] for x in columns]
            rows.append(row)
            
        line = map(lambda  x : '\t'.join(x), rows)
        return '\n'.join(lines)
        
    def get_attribute_as_quantity(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            raise AttributeError("unknown attribute '{0}'".format(attribute))
        
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        return attribute_values.values | attribute_values.unit
        
    
    def get_attributes_as_vector_quantities(self, attributes):
        for attribute in attributes:
            if not attribute in self.mapping_from_attribute_to_values_and_unit:
                raise AttributeError("unknown attribute '{0}'".format(attribute))
        
        unit_of_the_values = None
        results = []
        for attribute in attributes:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            if unit_of_the_values is None:
                unit_of_the_values = attribute_values.unit
                results.append(attribute_values.values)
            else:
                conversion_factor = attribute_values.unit.value_in(unit_of_the_values)
                if conversion_factor != 1.0:
                    results.append(attribute_values.values * conversion_factor)
                else:                    
                    results.append(attribute_values.values)
        
        results = numpy.dstack(results)[0]
        return results | unit_of_the_values
    
    
class Particles(object):
    class ScalarProperty(object):
        
        def  __init__(self, attribute_name):
            self.attribute_name = attribute_name
        
        def __get__(self, instance, owner):
            if instance == None:
                return self
            else:
                return instance.attributelist.get_attribute_as_quantity(self.attribute_name)
            
    class VectorProperty(object):
        
        def  __init__(self, attribute_names):
            self.attribute_names = attribute_names
        
        def __get__(self, instance, owner):
            if instance == None:
                return self
            else:
                return instance.attributelist.get_attributes_as_vector_quantities(self.attribute_names)
    
    """A set of particle objects"""
    def __init__(self, size = 0):
        self.particles = [Particle(i+1, self) for i in range(size)]
        self.attributelist = AttributeList(self.particles)
        self.history = []
        
    def __iter__(self):
        return self.attributelist.iter_particles()

    def __getitem__(self, index):
        return self.particles[index]
        
    def __len__(self):
        return len(self.particles)
        
    def savepoint(self):
        self.history.append(self.attributelist.copy())
    
    def get_timeline_of_attribute(self, particle, attribute):
        timeline = []
        for x in reversed(self.history):
            timeline.append((None, x.get_value_of(particle, attribute)))
        return timeline
            
    
            
        
class Stars(Particles):
    
    mass = Particles.ScalarProperty("mass")
    radius = Particles.ScalarProperty("radius")
    position = Particles.VectorProperty(["x","y","z"])
    velocity = Particles.VectorProperty(["vx","vy","vz"])
    
    
    def center_of_mass(self):
        masses, x_values, y_values, z_values = self.attributelist.get_values_of_all_particles_in_units(["mass","x","y","z"],[si.kg, si.m, si.m, si.m])
        total_mass = numpy.sum(masses)
        massx = numpy.sum(masses * x_values)
        massy = numpy.sum(masses * y_values)
        massz = numpy.sum(masses * z_values)
        position = numpy.array([massx, massy, massz])

        return values.new_quantity(position / total_mass, si.m)
    
    
    def center_of_mass_velocity(self):
        masses, x_values, y_values, z_values = self.attributelist.get_values_of_all_particles_in_units(["mass","vx","vy","vz"],[si.kg, si.m / si.s, si.m / si.s, si.m / si.s])
        total_mass = numpy.sum(masses)
        massx = numpy.sum(masses * x_values)
        massy = numpy.sum(masses * y_values)
        massz = numpy.sum(masses * z_values)
        position = numpy.array([massx, massy, massz])

        return values.new_quantity(position / total_mass, si.m / si.s)
        
        
        
            
        
            
class VectorAttribute(object):
    def __init__(self, names_of_the_scalar_components):
        self.names_of_the_scalar_components = names_of_the_scalar_components
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            return instance.set.attributelist.get_attribute_as_vector(instance, self.names_of_the_scalar_components)
    
    
    def __set__(self, instance, quantity):
        instance.set.attributelist.set_attribute_as_vector(instance, self.names_of_the_scalar_components, quantity)
            
class Particle(object):
    """A physical object or a physical region simulated as a 
    physical object (cload particle).
    
    All attributes defined on a particle are specific for 
    that particle (for example mass or position). A particle contains 
    a set of attributes, some attributes are *generic* and applicaple
    for multiple modules. Other attributes are *specific* and are 
    only applicable for a single module.
    """
    __slots__ = ["id", "attributes", "set", "_module_ids_to_index"]
    
    position = VectorAttribute(("x","y","z"))
    velocity = VectorAttribute(("vx","vy","vz"))
    
    def __init__(self, id = -1, set = None, **keyword_arguments):
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "attributes", {})
        object.__setattr__(self, "_module_ids_to_index", {})
        object.__setattr__(self, "set", set)
        
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if name_of_the_attribute == 'position':
            type(self).position.__set__(self, new_value_for_the_attribute)
            return
        if name_of_the_attribute == 'velocity':
            type(self).velocity.__set__(self, new_value_for_the_attribute)
            return
            
        if isinstance(new_value_for_the_attribute, values.Quantity):
            self.set.attributelist.set_value_of(self, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise Exception("attribute "+name_of_the_attribute+" does not have a valid value, values must have a unit")
    
    def __getattr__(self, name_of_the_attribute):
         return self.set.attributelist.get_value_of(self, name_of_the_attribute)
         
                
    def __str__(self):
        output = 'Particle '
        output += str(self.id)
        output += ''
        output += '\n'
        for name, value in self.set.attributelist.iter_values_of_particle(self):
            output += name
            output += ': {'
            output += str(value)
            output += '}, '
            output += '\n'
        return output
        
    def set_value_of_attribute(self, attribute, value, time = None):
        getattr(self, attribute).set_value_at_time(time, values[index])
        
    def set_default(self, attribute, quantity):
        if not attribute in self.set.attributelist.attributes():
            self.set.attributelist.set_value_of(self, attribute, quantity)
            
class Star(Particle):
    pass
