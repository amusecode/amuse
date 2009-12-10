"""
"""

from amuse.support.data import values
from amuse.support.units import si
from amuse.support.units import units

import numpy

class BasicUniqueKeyGenerator(object):
    
    def __init__(self):
        self.lowest_unique_key = 1
    
    def next(self):
        new_key = self.lowest_unique_key
        self.lowest_unique_key += 1
        return new_key
        
    def next_set_of_keys(self, length):
        if length == 0:
            return  []
            
        from_key = self.lowest_unique_key
        to_key = from_key + length;
        self.lowest_unique_key += length
        return numpy.arange(from_key, to_key)
        
UniqueKeyGenerator = BasicUniqueKeyGenerator()

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
    
    def __init__(self, particle_keys, attributes = [], lists_of_values = [], units = [], model_times = None):
        
        
        model_times = self._convert_model_times(model_times, len(particle_keys))
        
        self.mapping_from_attribute_to_values_and_unit = {}
        for attribute, values, unit in zip(attributes, lists_of_values, units):
            self.mapping_from_attribute_to_values_and_unit[attribute] = AttributeValues(
                attribute,
                unit,
                values,
                model_times
            )
        
        self.particle_keys = numpy.array(particle_keys)
        self.reindex()

    def _set_particles(self, keys, attributes = [], values = []):
        if len(values) != len(attributes):
            raise Exception(
                "you need to provide the same number of value list as attributes, found {0} attributes and {1} list of values".format(
                    len(attributes), len(values)
                )
            )
        if len(values) > 0 and len(keys) != len(values[0]):
            raise Exception(
                "you need to provide the same number of values as particles, found {0} values and {1} particles".format(
                    len(values[0]), len(keys)
                )
            )
        
        self.mapping_from_attribute_to_values_and_unit = {}
        
        for attribute, quantity in zip(attributes, values):
            self.mapping_from_attribute_to_values_and_unit[attribute] = AttributeValues(
                attribute,
                quantity.unit,
                quantity.number,
                None
            )
         
        self.particle_keys = numpy.array(keys)
        self.reindex()
    
    def _get_values(self, particles, attributes):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute in attributes:
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             if indices is None:
                 selected_values = attribute_values.values
             else:
                 selected_values = attribute_values.values.take(indices)
             results.append(attribute_values.unit.new_quantity(selected_values))
        
        return results
        
    def _set_values(self, particles, attributes, list_of_values_to_set, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, len(indices))
        
        previous_model_times = None
        results = []
        for attribute, values_to_set in zip(attributes, list_of_values_to_set):
            if attribute in self.mapping_from_attribute_to_values_and_unit:
                attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            else:
                attribute_values = AttributeValues(
                   attribute,
                   values_to_set.unit,
                   length = len(self.particle_keys)
                )
                self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
                 
            selected_values = values_to_set.value_in(attribute_values.unit)
            
            attribute_values.values.put(indices, selected_values)
            if not model_times is None:
                if not previous_model_times is attribute_values.model_times:
                    attribute_values.model_times.put(indices, model_times)
                    previous_model_times = attribute_values.model_times
            
        return results
    
    
    def _get_attributes(self):
        return sorted(self.mapping_from_attribute_to_values_and_unit.keys())
    
    
    def _has_key(self, key):
        return key in self.mapping_from_particle_to_index
        
    def _get_keys(self):
        return self.particle_keys
        
    def __len__(self):
        return len(self.particle_keys)
        
    def copy(self):
        copy = AttributeList([])
        copy.mapping_from_particle_to_index = self.mapping_from_particle_to_index.copy()
        copy.particle_keys = self.particle_keys.copy()
        for attribute, attribute_values in self.mapping_from_attribute_to_values_and_unit.iteritems():
            copy.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values.copy()
        return copy
        
    def get_value_of(self, particle_key, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            raise AttributeError("particle does not have a "+attribute)
        
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        index = self.mapping_from_particle_to_index[particle_key]
        
        return attribute_values.unit.new_quantity(attribute_values.values[index])
        
    def get_values_of_attributes_of_particle(self, particle_key, attributes):
        
        index = self.mapping_from_particle_to_index[particle_key]
        result = []
        for attribute in attribute:
            if not attribute in self.mapping_from_attribute_to_values_and_unit:
                raise AttributeError("particle does not have a "+attribute)
            
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            
            result.append(attribute_values.values[index] | attribute_values.unit)
        return result
        
    
    def set_value_of(self, particle_key, attribute, quantity):
        index = self.mapping_from_particle_to_index[particle_key]
            
        if index is None:
            raise Exception("particle with key '{0}' is not in this set".format(particle_key))
            
        attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
             
        value_to_set = quantity.value_in(attribute_values.unit)
        attribute_values.values[index] = value_to_set
        
            
    def iter_values_of_particle(self, particle_key):
        index = self.mapping_from_particle_to_index[particle_key]
        for attribute in self.mapping_from_attribute_to_values_and_unit:
            attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
            yield attribute, (attribute_values.values[index] | attribute_values.unit)
    
    
            
    def iter_values_of(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            return
            
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        values = attribute_values.values
        unit = attribute_values.unit
        particles = self.particle_keys
        
        for index in range(len(self.particle_keys)):
            yield particles[i], (values[i] | unit)
            
    def iter_particles(self):
        class ParticleView(object):
            __slots__=['key', 'index', 'version']
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
        for key in self.particle_keys:
            p = Particle(key, self)
            yield p
            index += 1
            
            
    
    def get_indices_of(self, particles):
        if particles is None:
            return numpy.arange(0,len(self.particle_keys))
            
        mapping_from_particle_to_index = self.mapping_from_particle_to_index 
        result = numpy.zeros(len(particles),dtype='int32')
        #result = [mapping_from_particle_to_index[id(particle)] for particle in particles]
        
        index = 0
        for index, particle_key in enumerate(particles):
            result[index] = mapping_from_particle_to_index[particle_key]
            index += 1
        return result
    
    def get_values_of_particles_in_units(self, particles, attributes, target_units):
        indices = self.get_indices_of(particles)
            
        results = []
        for attribute, target_unit in zip(attributes, target_units):
             attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
             value_of_unit_in_target_unit = attribute_values.unit.value_in(target_unit )
             if indices is None:
                 selected_values = attribute_values.values
             else:
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
    
    def _convert_model_times(self, model_times, length):
        if not model_times is None and isinstance(model_times, values.ScalarQuantity):
            return model_times.unit.new_quantity(numpy.linspace(model_times.number, model_times.number, length) )
        else:
            return model_times
    
    def set_values_of_particles_in_units(self, particles, attributes, list_of_values_to_set, source_units, model_times = None):
        indices = self.get_indices_of(particles)
        
        model_times = self._convert_model_times(model_times, len(indices))
        
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
                   length = len(self.particle_keys)
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
            
                
        other.set_values_of_particles_in_units(self.particle_keys, source_attributes, source_values, source_units)
        
    def remove_particles(self, particles):
        indices = self.get_indices_of(particles)
        
        mapping_from_attribute_to_values_and_unit = self.mapping_from_attribute_to_values_and_unit.copy()
        for attribute, attribute_values in mapping_from_attribute_to_values_and_unit.iteritems():
            attribute_values.values = numpy.delete(attribute_values.values,indices)
        
        self.particle_keys = numpy.delete(self.particle_keys,indices)
        self.reindex()
        
    def reindex(self):
        d = {}
        index = 0
        for particle_key in self.particle_keys:
            d[particle_key] = index
            index += 1
          
        self.mapping_from_particle_to_index = d
        
    def get_attribute_as_vector(self, particle_key, names_of_the_scalar_attributes):
        index = self.mapping_from_particle_to_index[particle_key]
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
                length = len(self.particle_keys)
            )
            self.mapping_from_attribute_to_values_and_unit[attribute] = attribute_values
        return attribute_values
        
    def set_attribute_as_vector(self, particle_key, names_of_the_scalar_attributes, quantity):
        index = self.mapping_from_particle_to_index[particle_key]
        vector = []
        for i, attribute  in enumerate(names_of_the_scalar_attributes):
            attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
        
            value_to_set = quantity[i].value_in(attribute_values.unit)
            attribute_values.values[index] = value_to_set
            
    def attributes(self):
        return set(self.mapping_from_attribute_to_values_and_unit.keys())
    
    def __str__(self):
        attributes = sorted(self.attributes())
        
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
        if len(self.particle_keys) > 40:
            values_to_show = list(self.particle_keys[:20])
            values_to_show.append(self.particle_keys[-20:])
        else:
            values_to_show = self.particle_keys
                    
        for value in values_to_show:
            column.append(str(value))
            
        column.append('========')
            
        rows = []
        for i in range(len(columns[0])):
            row = [x[i] for x in columns]
            rows.append(row)
            
        lines = map(lambda  x : '\t'.join(x), rows)
        return '\n'.join(lines)
        
    def get_attribute_as_quantity(self, attribute):
        if not attribute in self.mapping_from_attribute_to_values_and_unit:
            raise AttributeError("unknown attribute '{0}'".format(attribute))
        
        attribute_values = self.mapping_from_attribute_to_values_and_unit[attribute]
        
        return attribute_values.unit.new_quantity(attribute_values.values)
        
    def set_attribute_as_quantity(self, attribute, quantity):
        if not isinstance(quantity, values.VectorQuantity):
            raise AttributeError("can only set a vector of values")
        
        if not len(quantity) == len(self.particle_keys):
            raise AttributeError("vector of values must have the same length as the particles in the system")
            
        if not (
            isinstance(quantity._number[0], int) or isinstance(quantity._number[0], float)
            ):
            raise AttributeError("values must be ints or floats")
            
        attribute_values = self.get_attributevalues_for(attribute, quantity.unit)
        attribute_values.values = numpy.array(quantity._number)
        attribute_values.unit = quantity.unit
        
        
        
        
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
        return unit_of_the_values.new_quantity(results)
        
    def set_model_time(self, value): 
        model_times = self._convert_model_times(value, len(self.particle_keys))
        for attribute_values in self.mapping_from_attribute_to_values_and_unit.values():
            attribute_values.model_times = model_times
    
class Particles(object):
    class ScalarProperty(object):
        
        def  __init__(self, attribute_name):
            self.attribute_name = attribute_name
        
        def __get__(self, instance, owner):
            if instance == None:
                return self
            else:
                return instance._get_values(None, [self.attribute_name])[0]
        
        def __set__(self, instance, value):
            if instance == None:
                return self
            else:
                return instance._set_values(None, [self.attribute_name], [value])
                
    class VectorProperty(object):
        
        def  __init__(self, attribute_names):
            self.attribute_names = attribute_names
        
        def __get__(self, instance, owner):
            if instance == None:
                return self
            else:
                return instance.attributelist.get_attributes_as_vector_quantities(self.attribute_names)
                
        def __set__(self, instance, value):
            if instance == None:
                return self
            else:
                vectors = value.number
                split = numpy.hsplit(vectors,len(self.attribute_names))
                list_of_values = []
                for i in range(len(self.attribute_names)):
                    values = value.unit.new_quantity(split[i].reshape(len(split[i])))
                    list_of_values.append(values)
                    
                instance._set_values(None, self.attribute_names, list_of_values)
                
    class ModelTimeProperty(object):
        
        def __get__(self, instance, owner):
            if instance == None:
                return self
            else:
                raise Exception("TBD")
                
            
        def __set__(self, instance, value):
            if instance == None:
                return self
            else:
                instance.attributelist.set_model_time(value)
        
    
    """A set of particle objects"""
    def __init__(self, size = 0):
        particle_keys = UniqueKeyGenerator.next_set_of_keys(size)
        self.attributelist = AttributeList([])
        self._set_particles(particle_keys)
        
        self.previous = None
    
        
    def __iter__(self):
        index = 0
        for key in self.attributelist._get_keys():
            p = Particle(key, self)
            yield p
            index += 1

    def __getitem__(self, index):
        return Particle(self.attributelist._get_keys()[index], self)
        
    def __len__(self):
        return len(self.attributelist)
        
    def savepoint(self):
        instance = type(self)()
        instance.attributelist = self.attributelist.copy()
        instance.previous = self.previous
        self.previous = instance
        
    def iter_history(self):
        current = self
        while not current is None:
            yield current
            current = current.previous
    
    @property
    def history(self):
        return reversed(list(self.iter_history()))
        
    def get_timeline_of_attribute(self, particle_key, attribute):
        timeline = []
        for x in self.history:
            timeline.append((None, x.attributelist.get_value_of(particle_key, attribute)))
        return timeline

    model_time = ModelTimeProperty()
                    
    def copy(self):
         attributes = self.attributelist._state_attributes()
         keys = self._get_keys()
         values = self._get_values(None, attributes)
         result = Particles()
         result._set_particles(keys, attributes, values)
         return result
        
    def _set_particles(self, keys, attributes = [], values = []):
        self.attributelist._set_particles(keys, attributes, values)
    
    def _get_values(self, keys, attributes):
        return self.attributelist._get_values(keys, attributes)
        
    def _set_values(self, keys, attributes, values):
        self.attributelist._set_values(keys, attributes, values)
    
    def _get_attributes(self):
        return self.attributelist._get_attributes()
        
    def _get_keys(self):
        return self.attributelist._get_keys()
        
    def _has_key(self, key):
        return self.attributelist._has_key(key)
        
    def copy_values_of_attribute_to(self, attribute_name, particles):
        channel = self.new_channel_to(particles)
        channel.copy_attributes([attribute_name])
        
    def new_channel_to(self, other):
        return ParticleInformationChannel(self, other)
        
    def add_particles(self, particles):
        attributes = self._get_attributes()
        keys = particles._get_keys()
        values = particles._get_values(None, attributes)
        self._set_particles(keys, attributes, values)
        
    def copy_values_of_state_attributes_to(self, particles):
        channel = self.new_channel_to(particles)
        channel.copy_attributes(self.attributelist._state_attributes())  
              
class ParticleInformationChannel(object):
    
    def __init__(self, from_particles, to_particles):
        self.from_particles = from_particles
        self.to_particles = to_particles
        self._reindex()
        
    def _reindex(self):
        self.keys = self.intersecting_keys()
        #speed-up:
        #self.from_indices = self.from_particles._get_indices_of(self.keys)
        #self.to_indices = self.to_particles._get_indices_of(self.keys)
    
    def intersecting_keys(self):
        from_keys = self.from_particles._get_keys()
        return filter(lambda x : self.to_particles._has_key(x), from_keys)
        
    def copy_attributes(self, attributes):
        data = self.from_particles._get_values(self.keys, attributes)
        self.to_particles._set_values(self.keys, attributes, data)
    
class Stars(Particles):
    
    mass = Particles.ScalarProperty("mass")
    
    radius = Particles.ScalarProperty("radius")
    
    x = Particles.ScalarProperty("x")
    y = Particles.ScalarProperty("y")
    z = Particles.ScalarProperty("z")
    
    vx = Particles.ScalarProperty("vx")
    vy = Particles.ScalarProperty("vy")
    vz = Particles.ScalarProperty("vz")
    
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
        
    def kinetic_energy(self):
        mass = self.mass
        vx = self.vx
        vy = self.vy
        vz = self.vz
        v_squared = (vx * vx) + (vy * vy) + (vz * vz)
        m_v_squared = mass * v_squared
        return 0.5 * m_v_squared.sum()
        
    
    def potential_energy(self, smoothing_length_squared = 0 | units.m * units.m):
        mass = self.mass
        x_vector = self.x
        y_vector = self.y
        z_vector = self.z
        
        sum_of_energies = 0 | units.J
        
        for i in range(len(self)):
           x = x_vector[i]
           y = y_vector[i]
           z = z_vector[i]
           dx = x - x_vector[i+1:]
           dy = y - y_vector[i+1:]
           dz = z - z_vector[i+1:]
           dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
           dr = (dr_squared+smoothing_length_squared).sqrt()
           m_m = mass[i] * mass[i+1:]
           
           energy_of_this_particle = (m_m / dr).sum()
           sum_of_energies +=  -1 * units.G * energy_of_this_particle
            
        return sum_of_energies
        
        
        
        
        
            
        
            
class VectorAttribute(object):
    def __init__(self, names_of_the_scalar_components):
        self.names_of_the_scalar_components = names_of_the_scalar_components
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            return instance.particles_set.attributelist.get_attribute_as_vector(instance.key, self.names_of_the_scalar_components)
    
    
    def __set__(self, instance, quantity):
        instance.particles_set.attributelist.set_attribute_as_vector(instance.key, self.names_of_the_scalar_components, quantity)
            
class Particle(object):
    """A physical object or a physical region simulated as a 
    physical object (cload particle).
    
    All attributes defined on a particle are specific for 
    that particle (for example mass or position). A particle contains 
    a set of attributes, some attributes are *generic* and applicaple
    for multiple modules. Other attributes are *specific* and are 
    only applicable for a single module.
    """
    __slots__ = ["key", "particles_set"]
    
    
    position = VectorAttribute(("x","y","z"))
    velocity = VectorAttribute(("vx","vy","vz"))
    
    def __init__(self, key, particles_set = None, **keyword_arguments):
        object.__setattr__(self, "key", key)
        object.__setattr__(self, "particles_set", particles_set)
        
        for attribute_name in keyword_arguments:
            attribute_value = keyword_arguments[attribute_name]
            setattr(self, attribute_name, attribute_value)
            
    def __setattr__(self, name_of_the_attribute, new_value_for_the_attribute):
        if hasattr(type(self), name_of_the_attribute):
            getattr(type(self), name_of_the_attribute).__set__(self, new_value_for_the_attribute)
            return
            
        if isinstance(new_value_for_the_attribute, values.Quantity):
            self.particles_set.attributelist.set_value_of(self.key, name_of_the_attribute, new_value_for_the_attribute)
        else:
            raise Exception("attribute "+name_of_the_attribute+" does not have a valid value, values must have a unit")
    
    def __getattr__(self, name_of_the_attribute):
         return self.particles_set.attributelist._get_values([self.key], [name_of_the_attribute])[0][0]
         
                
    def __str__(self):
        output = 'Particle '
        output += str(self.key)
        output += ''
        output += '\n'
        for name, value in self.particles_set.attributelist.iter_values_of_particle(self.key):
            output += name
            output += ': {'
            output += str(value)
            output += '}, '
            output += '\n'
        return output
        
    def set_default(self, attribute, quantity):
        if not attribute in self.set.attributelist.attributes():
            self.particles_set.attributelist.set_value_of(self, attribute, quantity)
            
    def get_timeline_of_attribute(self, attribute):
        return self.particles_set.get_timeline_of_attribute(self.key, attribute)
            
