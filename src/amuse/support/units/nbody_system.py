"""
"""
from amuse.support.units import units
from amuse.support.units import core

import numpy

class nbody_unit(core.unit):
    def __init__(self, unit):
        self.unit = unit
    def __str__(self):
        return 'nbody '+self.unit.quantity
    @property
    def base(self):
        return ((1,self),)
    @property
    def factor(self):
        return 1.0

length = nbody_unit(units.m)
time = nbody_unit(units.s)
mass = nbody_unit(units.kg)
 
speed = length / time

class nbody_to_si(object): 
    
    def __init__(self, value1 , value2):
        self.value1 = value1
        self.value2 = value2
        if self.unit1.base == self.unit2.base:
            raise Exception("Must provide two orthogonal units for example mass and length or time and length")
        if len(self.unit1.base) > 1:
            raise Exception("Currently cannot handle unit with more than one base, only meters or seconds or kg")
        if len(self.unit2.base) > 1:
            raise Exception("Currently cannot handle unit with more than one base, only meters or seconds or kg")
        self.Gis1 = units.one / units.G
        
    @property
    def unit1(self):
        return self.value1.unit 
        
    @property
    def unit2(self):
        return self.value2.unit 
    
    @property    
    def units(self):
        vector = [0.0,0.0,0.0]
        unit_to_index = {}
        for i, x in enumerate(self.Gis1.unit.base):
            n , unit = x
            unit_to_index[unit] = i
    
        exponents_of_the_bases =  numpy.zeros((3,3))
        factors_of_the_bases =  numpy.zeros(3)
        for row, value in enumerate((self.Gis1, self.value1, self.value2)):
            for n, unit in value.unit.base:
                exponents_of_the_bases[row, unit_to_index[unit]] = n
            factors_of_the_bases[row] = value.number * value.unit.factor  
        log_factors_of_the_bases = numpy.log(factors_of_the_bases)
        conversion_factors = numpy.exp(numpy.linalg.solve(exponents_of_the_bases, log_factors_of_the_bases))
        
        result = []
        nbody_units = mass, length, time
        
        for n , unit  in self.Gis1.unit.base:
            index = unit_to_index[unit]
            conversion_factor_for_this_base_unit = conversion_factors[index]
            for nbody_unit in nbody_units:
                if nbody_unit.unit == unit:
                    result.append((nbody_unit, conversion_factor_for_this_base_unit * unit))
        return result
    def to_si(self, value):
        nbody_units_in_si = self.units
        base = value.unit.base
        factor = value.unit.factor
        number = value.number
        new_unit = 1
        for n, unit in base:
            found = False
            for unit_nbody,unit_in_si in nbody_units_in_si:
                if unit_nbody==unit:
                    factor = factor * (unit_in_si.factor ** n)
                    new_unit = new_unit * (unit_in_si.base[0][1] ** n)
                    found = True
                    break
            if not found:
                new_unit = new_unit * (unit ** n)
        return new_unit(number * factor)
        
        
        
    def to_nbody(self, value):
        nbody_units_in_si = self.units
        base = value.unit.base
        factor = value.unit.factor
        number = value.number
        new_unit = 1
        for n, unit in base:
            found = False
            for unit_in_nbody, unit_in_si in nbody_units_in_si:
                base_in_si = unit_in_si.base[0][1]
                if base_in_si == unit:
                    factor = factor / (unit_in_si.factor ** n)
                    new_unit = new_unit * (unit_in_nbody.base[0][1] ** n)
                    found = True
                    break;
            if not found:
                new_unit = new_unit * (unit ** n)
        return new_unit(number * factor)
        
                
                    
        

       
