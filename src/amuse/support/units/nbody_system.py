from amuse.support.units import units
from amuse.support.units import core
from amuse.support.units import constants
from amuse.support.units import generic_unit_converter
from amuse.support import exceptions

from amuse.support.data.values import new_quantity


import numpy

"""
"""
class nbody_unit(core.base_unit):
    def __init__(self, unit_in_si, system):
        core.base_unit.__init__(self, unit_in_si.quantity, unit_in_si.name, unit_in_si.symbol, system)
        self.unit_in_si = unit_in_si
        
    def __str__(self):
        return 'nbody '+self.unit_in_si.quantity
        
  

    def is_generic(self):
        return True
    
    
nbody_system = core.system('nbody')

length = nbody_unit(units.m, nbody_system)
time = nbody_unit(units.s, nbody_system)
mass = nbody_unit(units.kg, nbody_system)
acceleration = length / (time ** 2)
potential = (length ** 2) / (time ** 2)
energy = mass * potential
specific_energy = potential
speed = length / time
volume = (length ** 3)
density = mass / volume
momentum_density = density * speed
energy_density = density * specific_energy
G = 1 | (length**3) / (mass * (time**2))

def is_nbody_unit(unit):
    for factor, x in unit.base:
        if x.is_generic():
            return True
    return False
    
    
class nbody_to_si(generic_unit_converter.ConvertBetweenGenericAndSiUnits): 
    def __init__(self, value1 , value2):
        generic_unit_converter.ConvertBetweenGenericAndSiUnits.__init__(self,constants.G, value1, value2)
        
    def to_nbody(self, value):
        return self.to_generic(value)
        
    
    def _old_unit_to_unit_in_nbody(self, unit):
        nbody_units_in_si = self.units
        base = unit.base
        factor = unit.factor
        new_unit = 1
        for n, unit in base:
            unit_in_nbody, unit_in_si = self.find_nbody_unit_for(unit)
            if not unit_in_si is None:
                factor = factor / (unit_in_si.factor ** n)
                new_unit *= (unit_in_nbody.base[0][1] ** n)
            else:
                new_unit *= (unit ** n)
        return factor * new_unit
        
    def as_converter_from_si_to_nbody(self):
        class SiToNBodyConverter(object):
            def __init__(self, nbody_to_si):
                self.nbody_to_si = nbody_to_si
            
            def from_source_to_target(self, quantity):
                if hasattr(quantity, 'unit'):
                    return self.nbody_to_si.to_nbody(quantity) 
                else:
                    return quantity
                
            def from_target_to_source(self, quantity):
                if hasattr(quantity, 'unit'):
                    return self.nbody_to_si.to_si(quantity)
                else:
                    return quantity
                    
                
        return SiToNBodyConverter(self)
        
    def as_converter_from_nbody_to_si(self):
        
        class NBodyToSiConverter(object):
            def __init__(self, nbody_to_si):
                self.nbody_to_si = nbody_to_si
            
            def from_source_to_target(self, quantity):
                if hasattr(quantity, 'unit'):
                    return self.nbody_to_si.to_si(quantity) 
                else:
                    return quantity
                
            def from_target_to_source(self, quantity):
                if hasattr(quantity, 'unit'):
                    return self.nbody_to_si.to_nbody(quantity)
                else:
                    return quantity
                    
                
        return NBodyToSiConverter(self)


    def as_converter_from_si_to_generic(self):
        class SiToNBodyConverter(object):
            def __init__(self, nbody_to_si):
                self.nbody_to_si = nbody_to_si
            
            def from_source_to_target(self, quantity):
                if hasattr(quantity, 'unit'):
                    return self.nbody_to_si.to_nbody(quantity)
                else:
                    return quantity
                
            def from_target_to_source(self, quantity):
                if hasattr(quantity, 'unit'):
                    return self.nbody_to_si.to_si(quantity)
                else:
                    return quantity
                    
                
        return SiToNBodyConverter(self)
        
    @property
    def units(self):
        conversion_factors = self.conversion_factors()
        result = []
        generic_units = mass, length, time #, temperature, current, luminous_intensity

        for n, unit  in enumerate(self.list_of_available_units):
            conversion_factor_for_this_base_unit = conversion_factors[n]
            for generic_unit in generic_units:
                if generic_unit.unit_in_si == unit:
                    result.append((generic_unit, conversion_factor_for_this_base_unit * unit))

        return result
