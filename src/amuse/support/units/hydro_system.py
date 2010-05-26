
"""
"""
from amuse.support.units import units
from amuse.support.units import core
from amuse.support.units import constants

from amuse.support.data.values import new_quantity

import numpy

class hydro_unit(core.base_unit):
    def __init__(self, unit_in_si, system):
        core.base_unit.__init__(self, unit_in_si.quantity, unit_in_si.name, unit_in_si.symbol, system)
        self.unit_in_si = unit_in_si

    def __str__(self):
        return 'hydro '+self.unit_in_si.quantity

hydro_system = core.system('hydro')

length = hydro_unit(units.m, hydro_system)
time = hydro_unit(units.s, hydro_system)
mass = hydro_unit(units.kg, hydro_system)
acceleration = length / (time ** 2)
force = mass*acceleration
potential = (length ** 2) / (time ** 2)
energy = mass * potential
speed = length / time

def is_hydro_unit(unit):
    for factor, x in unit.base:
        if isinstance(x, hydro_unit):
            return True
    return False

class hydro_to_si(object):
    DEFAULT_CONVERTER = None

    def __init__(self, value1 , value2, value3):
        self.value1 = value1
        self.value2 = value2
        self.value3 = value3
        
        self.new_base = numpy.mat(numpy.zeros((3,3)))
        self.new_base_inv = numpy.mat(numpy.zeros((3,3)))

        #isolate the 3 base units
        available_units = set()
        for unit in [value1.unit, value2.unit, value3.unit]:
            for i in unit.base:
                available_units.add(i[1])
        if not len(available_units) is 3:
            raise Exception("Must provide three orthogonal units")
        self.list_of_available_units = list(available_units)

        self.new_base = self.determine_new_base()
        if self.matrixrank(self.new_base) < 3:
            raise Exception("Must provide three orthogonal units, e.g. mass,length,time or mass, velocity, time")
        self.new_base_inv = self.new_base ** -1

        self.set_default_converter_if_uninitialised(self)

    def matrixrank(self, A, tol=1e-8):
        s = numpy.linalg.svd(A,compute_uv=0)
        return numpy.sum(numpy.where( s>tol, 1, 0 ) )

    @property
    def unit1(self):
        return self.value1.unit

    @property
    def unit2(self):
        return self.value2.unit

    @property
    def unit3(self):
        return self.value3.unit

    def determine_new_base(self):
        matrix = numpy.mat(numpy.zeros((3,3)))
        for row, value in enumerate((self.value1, self.value2, self.value3)):
            for n, unit in value.unit.base:
                matrix[row, [i for i, j in enumerate(self.list_of_available_units) if j == unit]] = n
        return matrix

    def conversion_factors(self):
        factors_of_the_bases =  numpy.mat(numpy.zeros((3,1)))
        for row, value in enumerate((self.value1, self.value2, self.value3)):
            factors_of_the_bases[row] = value.number * value.unit.factor

        log_factors_of_the_bases = numpy.log(factors_of_the_bases)
        return numpy.exp(self.new_base_inv*log_factors_of_the_bases)

    @property
    def units(self):
        conversion_factors = self.conversion_factors()
        result = []
        hydro_units = mass, length, time

        for n, unit  in enumerate(self.list_of_available_units):
            conversion_factor_for_this_base_unit = conversion_factors[n]
            for hydro_unit in hydro_units:
                if hydro_unit.unit_in_si == unit:
                    result.append((hydro_unit, conversion_factor_for_this_base_unit * unit))

        return result

    def find_si_unit_for(self, unit):
        for unit_hydro, unit_in_si in self.units:
            if unit_hydro == unit:
                return unit_hydro, unit_in_si
        return None, None

    def find_hydro_unit_for(self, unit):
        for unit_hydro, unit_in_si in self.units:
            base_in_si = unit_in_si.base[0][1]
            if base_in_si == unit:
                return unit_hydro, unit_in_si
        return None, None

    def to_si(self, value):
        factor = value.unit.factor
        number = value.number
        new_unit = 1
        base = value.unit.base
        
        if not base:
            return value
        
        for n, unit in base:
            unit_in_hydro, unit_in_si = self.find_si_unit_for(unit)
            if not unit_in_si is None:
                factor = factor * (unit_in_si.factor ** n)
                new_unit *= (unit_in_si.base[0][1] ** n)
            else:
                new_unit *= (unit ** n)
        return new_quantity(number * factor, new_unit)

    def to_hydro(self, value):
        hydro_units_in_si = self.units
        base = value.unit.base
        factor = value.unit.factor
        number = value.number
        new_unit = 1
        
        if not base:
            return value
            
        for n, unit in base:
            unit_in_hydro, unit_in_si = self.find_hydro_unit_for(unit)
            if not unit_in_si is None:
                factor = factor / (unit_in_si.factor ** n)
                new_unit *= (unit_in_hydro.base[0][1] ** n)
            else:
                new_unit *= (unit ** n)
        return new_quantity(number * factor, new_unit)

    def set_as_default(self):
        """Install this converter as the default converter for the
        modules. When a native nbody module is created it will choose
        this converter (if no other converter is given during creation)
        """
        self.set_default_converter(self)

    @classmethod
    def get_default(cls):
        if cls.DEFAULT_CONVERTER is None:
            raise Exception("Asked for the default nbody to SI converter,"
            " but no converter has been set!.\n"
            "Please create a nbody_to_si converter first,"
            " and use the 'set_as_default' method.")
        else:
            return cls.DEFAULT_CONVERTER

    @classmethod
    def set_default_converter(cls, object):
        cls.DEFAULT_CONVERTER = object

    @classmethod
    def set_default_converter_if_uninitialised(cls, object):
        if cls.DEFAULT_CONVERTER is None:
            cls.set_default_converter(object)

class noconvert_hydro_to_si(object):

    def __init__(self):
        pass

    def to_si(self, value):
        return value

    def to_hydro(self, value):
        return value
