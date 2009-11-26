from amuse.support.units import si
from amuse.support.units import units
from amuse.support.units import nbody_system

import numpy


class AttributeDefinition(object):
    def __init__(self, 
            unit, 
            name = None,
            names = None, 
            setup_parameters = None,
            state_parameters = None,
            getter=None, 
            setter=None, 
            description="",  
            default=None,
        ):
        
        self.name = name
        if not name is None and names is None:
            self.names = [name]
        else:
            self.names = names
            
        self.description = description
        self.unit = unit
        self.default_value = default
        self.getter = getter
        self.setter = setter
        self.setup_parameters = setup_parameters
        if state_parameters is None:
            self.state_parameters = self.setup_parameters
        else:
            self.state_parameters = state_parameters
            
    def is_required_for_setup(self):
        return not self.setup_parameters is None
        
        
    def for_setup_fill_arguments_for_attributelist_get(self, attributes, units, keywords):
        for name, parameter_name in zip(self.names, self.setup_parameters):
            attributes.append(name)
            units.append(self.unit)
            keywords.append(parameter_name)
        
    
    def for_state_fill_arguments_for_attributelist_get(self, attributes, units, keywords):
        for name, parameter_name in zip(self.names, self.state_parameters):
            attributes.append(self.name)
            units.append(self.unit)
            keywords.append(parameter_name)
            
    
    def for_state_fill_arguments_for_attributelist_set(self, states,  attributes, units, values):
        for name, parameter_name in zip(self.names, self.state_parameters):
            units.append(self.unit)
            attributes.append(name)
            values.append(states[parameter_name])
            
    
    def for_setter_fill_arguments_for_attributelist_get(self, attributes, units, keywords):
        for name, parameter_name in zip(self.names, self.setter[1]):
            attributes.append(self.name)
            units.append(self.unit)
            keywords.append(parameter_name)
            
    
    def for_getter_fill_arguments_for_attributelist_set(self, states,  attributes, units, values):
        for name, parameter_name in zip(self.names, self.getter[1]):
            units.append(self.unit)
            attributes.append(name)
            values.append(states[parameter_name])
         
        
        
    
        

class AttributeDefinition_Old(object):
    def __init__(self, name, description, unit, default_value):
        self.name = name
        self.description = description
        self.unit = unit
        self.default_value = default_value


class DomainMetaclass(type):
    def __new__(metaclass, name, bases, dict):
        replacement_dictionary = {}
        for key, value in dict.iteritems():
            if isinstance(value, tuple):
                default_value, description = value
                replacement_dictionary[key] = AttributeDefinition_Old(
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
