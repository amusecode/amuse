from amuse.units import generic_unit_system
from amuse.units.quantities import new_quantity

class ScalingConverter(object):

    def __init__(
        self,
        length = 1,
        time = 1,
        mass = 1,
        current = 1,
        temperature = 1,
        amount = 1,
        luminous_intensity = 1
    ):
        self.factors = {}
        
        self.factors[generic_unit_system.mass] = mass
        self.factors[generic_unit_system.length] = length
        self.factors[generic_unit_system.time] = time
        self.factors[generic_unit_system.temperature] = temperature
        self.factors[generic_unit_system.current] = current
        self.factors[generic_unit_system.luminous_intensity] = luminous_intensity
        
    def reversed(self):
        return ScalingConverter(
            length = 1 / self.factors[generic_unit_system.length],
            time = 1 / self.factors[generic_unit_system.time],
            mass =  1 / self.factors[generic_unit_system.mass],
            current = 1 / self.factors[generic_unit_system.current],
            temperature =  1 / self.factors[generic_unit_system.temperature],
            amount = 1,
            luminous_intensity = 1 / self.factors[generic_unit_system.luminous_intensity]
        )
        
    def convert(self, quantity):
        unit = quantity.unit
        value = quantity.value_in(unit)
        
        base = unit.base
        if not base:
            return quantity
            
        new_unit = 1
        factor = unit.factor
        
        for n, unit in base:
            if unit in self.factors:
                factor_for_unit = self.factors[unit]
                factor = factor * (factor_for_unit ** n)
                new_unit *= (unit ** n)
            else:
                new_unit *= (unit ** n)
        return new_quantity(value * factor, new_unit)