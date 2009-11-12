from amuse.support.units import si
from amuse.support.units import units
from amuse.support.units import nbody_system

import numpy

class AttributeDefinition(object):
    def __init__(self, name, description, unit, default_value):
        self.name = name
        self.description = description
        self.unit = unit
        self.default_value = default_value
    
    def convert_to_numbers_in_legacy_units(self, object, values):
        if nbody_system.is_nbody_unit(self.unit):
            values = [object.convert_nbody.to_nbody(x) for x in values]
        return [x.value_in(self.unit) for x in values]
    
    def convert_to_quantities(self, object, numbers):
        values = [self.unit(x) for x in numbers]
        if nbody_system.is_nbody_unit(self.unit):
            values = [object.convert_nbody.to_si(x) for x in values]
        return values
        
    def set_values(self, object, ids, values):
        numbers = self.convert_to_numbers_in_legacy_units(object, values)
        self.set_legacy_values(object, ids, numbers)
        
    def get_values(self, object, ids):
        numbers = self.get_legacy_values(object, ids)
        return self.convert_to_quantities(object, numbers)
        
class ScalarAttributeDefinition(AttributeDefinition):
    def __init__(self, set_method, get_method, parameter_name,  name, description, unit, default_value):
        AttributeDefinition.__init__(self, name, description, unit, default_value)
        self.set_method = set_method
        self.get_method = get_method
        self.parameter_name = parameter_name
        
    def set_keyword_arguments(self, object, values, keyword_arguments):
        numbers = self.convert_to_numbers_in_legacy_units(object, values)
        keyword_arguments[self.parameter_name] = numbers
    
    def get_keyword_results(self, object, keyword_results):
        numbers = keyword_results[self.parameter_name]
        values = self.convert_to_quantities(object, numbers)
        return values
        

class AttributeException(AttributeError):
    template = ("Could not {0} value(s) for atttribute '{1}' of a '{2}' object, got errorcode <{3}>")
        
    def __init__(self, object, parameter_name, errorcode, is_get):
        print is_get, "get" if is_get  else "set"
        AttributeError.__init__(self, self.template.format(
            "get" if is_get else "set",
            parameter_name,
            type(object).__name__,
            errorcode
        ))
        self.errorcode = errorcode
        self.parameter_name = parameter_name


class ScalarAttributeDefinition_Next(AttributeDefinition):
    def __init__(self, set_method, get_method, parameter_name,  name, description, unit, default_value):
        AttributeDefinition.__init__(self, name, description, unit, default_value)
        self.set_method = set_method
        self.get_method = get_method
        self.parameter_name = parameter_name
        
    def set_keyword_arguments(self, object, values, keyword_arguments):
        numbers = self.convert_to_numbers_in_legacy_units(object, values)
        keyword_arguments[self.parameter_name] = numbers
    
    def get_keyword_results(self, object, keyword_results):
        numbers = keyword_results[self.parameter_name]
        values = self.convert_to_quantities(object, numbers)
        return values
    
    def set_legacy_values(self, object, ids, values):
        error = getattr(object, self.set_method)(ids, values)
        if error < 0:
            raise AttributeException(object, self.name, error, False)
    
    def get_legacy_values(self, object, ids):
        values, error = getattr(object, self.get_method)(ids)
        if error < 0:
            raise AttributeException(object, self.name, error, True)
        return values
        
class VectorAttributeDefinition(AttributeDefinition):
    def __init__(self, set_method, get_method, parameter_names,  name, description, unit, default_value):
        AttributeDefinition.__init__(self, name, description, unit, default_value)
        self.set_method = set_method
        self.get_method = get_method
        self.parameter_names = parameter_names
        
    def set_keyword_arguments(self, object, values, keyword_arguments):
        numbers = self.convert_to_numbers_in_legacy_units(object, values)
        for i, x in enumerate(self.parameter_names):
            keyword_arguments[x] = [number[i] for number in numbers]
    
    
    def get_keyword_results(self, object, keyword_results):
        numbers = [None for x in self.parameter_names]
        for i, x in enumerate(self.parameter_names):
            numbers[i] = keyword_results[x]
        numbers = numpy.array(numbers)
        numbers = numbers.transpose()
        values = self.convert_to_quantities(object, numbers)
        return values

class DomainMetaclass(type):
    def __new__(metaclass, name, bases, dict):
        replacement_dictionary = {}
        for key, value in dict.iteritems():
            if isinstance(value, tuple):
                default_value, description = value
                replacement_dictionary[key] = AttributeDefinition(
                        key, description, 
                        default_value.unit, default_value)
            else:
                replacement_dictionary[key] = value
        return type.__new__(metaclass, name, bases, dict)
        
        
class Domain(object):
    __metaclass__ = DomainMetaclass
    time = 0.0 | si.s , "model time"
    
class Gravity(Domain):
    mass = 0.0 | si.kg , "the mass of a star"
    position = [0.0, 0.0, 0.0] | si.m , "the position vector of a star"
    velocity = [0.0, 0.0, 0.0] | si.m / si.s , "the velocity vector of a star"
    radius = 0.0 | si.m , "the radius of a star"
    acceleration = [0.0, 0.0, 0.0] | si.m / (si.s ** 2), "the acceleraration vector of a star"

class Hydrodynamics(Domain):
    pressure = 0.0 | units.Pa , "the pressure in a region of space"
    density = 0.0 | si.kg / (si.m ** 3), "the density of molecules or solid matter"
    temperature = 0.0 | si.K , "the temperature of the gas"
    magnetic_field = 0.0 | units.tesla, "magnetic field created by gas and stars"
    velovity_field = 0.0 | si.m / si.s  , "velocity of the gas"
    gravity_potential = 0.0 | si.no_unit  , "gravity forces from stars and gas"
    viscosity = 0.0 | si.no_unit, "viscosity of the gas cloud"
    
class RadiativeTransfer(Domain):
    temperature_gas = 0.0 | si.K , "the temperature of the gas"
    temperature_dust = 0.0 | si.K , "the temperature of the dust"
    temperature_background = 0.0 | si.K , "the temperature of the background"
    density = 0.0 | si.mol / (si.m**3), "modulecular density"
    magnetic_field = 0.0 | units.tesla, "magnetic field created by gas and stars"
    velovity_field = 0.0 | si.m / si.s  , "velocity of the gas"
    
    
    
class StellarEvolution(Domain):
    mass = 0.0 | si.kg , "the mass of a star"
    radius = 0.0 | si.m , "the radius of a star"
    age = 0.0 | si.s , "the age of a star, time evolved since star formation"


class SseCode(StellarEvolution):
    zams_mass = 0.0 | si.kg , "the mass of a star after formation"
    type = 0 | si.no_unit, "stars evolve through typical stages, during each stage one can classify a star as belonging to a specific type"
    luminosity = 0.0 | si.cd / (si.m ** 2), "brightness of a star"
    radius = 0.0 | si.m, "total radius of a star"
    core_mass = 0.0 | si.kg, "mass of the innermost layer of a star"
    core_radius = 0.0 | si.m, "radius of the innermost layer of a star"
    envelope_mass = 0.0 | si.kg, "mass of the radiative / convective envelope around the core of a star"
    envelope_radius = 0.0 | si.m, "radius of the radiative / convective envelope around the core of a star"
    spin = 0.0 | si.m / si.s, "speed of rotation around the central axis of a star"
    epoch = 0.0 | si.s, "set when a star changes type"
    physical_time = 0.0 | si.s, "age of a star relative to last change of type"
