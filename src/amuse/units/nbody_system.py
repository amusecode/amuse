"""

The n-body unit system knows the three base quantities in the
International System of Quantities, I.S.Q. and defines 
the gravitational constant to be 1:

G = 1 | (length**3) / (mass * (time**2))

+-------------------+-----------------------------------+-----------------+
|Base quantity      |Name in generic unit               |Name in S.I. unit|
+-------------------+-----------------------------------+-----------------+
|length             |generic_system.length              |units.m          |
+-------------------+-----------------------------------+-----------------+
|time               |generic_system.time                |units.s          |
+-------------------+-----------------------------------+-----------------+
|mass               |generic_system.mass                |units.kg         |
+-------------------+-----------------------------------+-----------------+

Derived quantities
~~~~~~~~~~~~~~~~~~

+------------------+--------------------------------+
|acceleration      |length / (time ** 2)            |
+------------------+--------------------------------+
|potential         |(length ** 2) / (time ** 2)     |
+------------------+--------------------------------+
|energy            |mass * potential                |
+------------------+--------------------------------+
|specific_energy   |potential                       |
+------------------+--------------------------------+
|speed             |length / time                   |
+------------------+--------------------------------+
|volume            |(length ** 3)                   |
+------------------+--------------------------------+
|density           |mass / volume                   |
+------------------+--------------------------------+
|momentum_density  |density * speed                 |
+------------------+--------------------------------+
|energy_density    |density * specific_energy       |
+------------------+--------------------------------+


"""
from amuse.units import units
from amuse.units import core
from amuse.units import constants
from amuse.units import generic_unit_converter
from amuse.units import generic_unit_system
from amuse.units.quantities import new_quantity
from amuse.support import exceptions

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
    
length = generic_unit_system.length
time =  generic_unit_system.time
mass =  generic_unit_system.mass

acceleration = length / (time ** 2)
potential = (length ** 2) / (time ** 2)
energy = mass * potential
specific_energy = potential
speed = length / time
volume = (length ** 3)
density = mass / volume
pressure = mass / length / (time ** 2)
momentum_density = density * speed
energy_density = density * specific_energy
G = 1. | (length**3) / (mass * (time**2))

def is_nbody_unit(unit):
    for factor, x in unit.base:
        if x.is_generic():
            return True
    return False
    
class SiToNBodyConverter(object):
    def __init__(self, nbody_to_si):
        self.nbody_to_si = nbody_to_si
    
    def from_source_to_target(self, quantity):
        if hasattr(quantity, 'unit') and not quantity.unit.is_non_numeric():
            return self.nbody_to_si.to_nbody(quantity)
        else:
            return quantity
        
    def from_target_to_source(self, quantity):
        if hasattr(quantity, 'unit') and not quantity.unit.is_non_numeric():
            return self.nbody_to_si.to_si(quantity)
        else:
            return quantity
    
class nbody_to_si(generic_unit_converter.ConvertBetweenGenericAndSiUnits): 
    def __init__(self, value1, value2):
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
        return self.as_converter_from_si_to_generic()
        
    def as_converter_from_nbody_to_si(self):
        return self.as_converter_from_generic_to_si()


    def as_converter_from_si_to_generic(self):
        return SiToNBodyConverter(self)
        
