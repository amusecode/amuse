import numpy
from amuse.support.units.generic_unit_system import *

class generic_to_si(object):
    DEFAULT_CONVERTER = None

    def __init__(self, *kargs):
        self.values = kargs
        self.units_of_values = [i.unit for i in self.values]
        self.system_rank = len(kargs)
        
        self.new_base = numpy.mat(numpy.zeros((self.system_rank,self.system_rank)))
        self.new_base_inv = numpy.mat(numpy.zeros((self.system_rank,self.system_rank)))

        available_units = set()
        for unit in self.units_of_values:#[self.value.unit, self.value2.unit, self.value3.unit]:
            for i in unit.base:
                available_units.add(i[1])
        if not len(available_units) is self.system_rank:
            raise Exception("Must provide three orthogonal units")
        self.list_of_available_units = list(available_units)

        self.new_base = self.determine_new_base()
        if self.matrixrank(self.new_base) < self.system_rank:
            raise Exception("Must provide three orthogonal units, e.g. mass,length,time or mass, velocity, time")
        self.new_base_inv = self.new_base ** -1
        print self.new_base
        print self.new_base_inv
        self.set_default_converter_if_uninitialised(self)

    def matrixrank(self, A, tol=1e-8):
        s = numpy.linalg.svd(A,compute_uv=0)
        return numpy.sum(numpy.where( s>tol, 1, 0 ) )

    def determine_new_base(self):
        matrix = numpy.mat(numpy.zeros((self.system_rank,self.system_rank)))
        for row, value in enumerate(self.values):
            for n, unit in value.unit.base:
                matrix[row, [i for i, j in enumerate(self.list_of_available_units) if j == unit]] = n
        return matrix

    def conversion_factors(self):
        factors_of_the_bases =  numpy.mat(numpy.zeros((self.system_rank,1)))
        for row, value in enumerate(self.values):
            factors_of_the_bases[row] = value.number * value.unit.factor

        log_factors_of_the_bases = numpy.log(factors_of_the_bases)
        return numpy.exp(self.new_base_inv*log_factors_of_the_bases)

    @property
    def units(self):
        conversion_factors = self.conversion_factors()
        result = []
        generic_units = mass, length, time, temperature, current, luminous_intensity

        for n, unit  in enumerate(self.list_of_available_units):
            conversion_factor_for_this_base_unit = conversion_factors[n]
            for generic_unit in generic_units:
                if generic_unit.unit_in_si == unit:
                    result.append((generic_unit, conversion_factor_for_this_base_unit * unit))

        return result

    def find_si_unit_for(self, unit):
        for unit_generic, unit_in_si in self.units:
            if unit_generic == unit:
                return unit_generic, unit_in_si
        return None, None

    def find_generic_unit_for(self, unit):
        for unit_generic, unit_in_si in self.units:
            base_in_si = unit_in_si.base[0][1]
            if base_in_si == unit:
                return unit_generic, unit_in_si
        return None, None

    def to_si(self, value):
        factor = value.unit.factor
        number = value.number
        new_unit = 1
        base = value.unit.base
        
        if not base:
            return value
        
        for n, unit in base:
            unit_in_generic, unit_in_si = self.find_si_unit_for(unit)
            if not unit_in_si is None:
                factor = factor * (unit_in_si.factor ** n)
                new_unit *= (unit_in_si.base[0][1] ** n)
            else:
                new_unit *= (unit ** n)
        return new_quantity(number * factor, new_unit)

    def to_generic(self, value):
        generic_units_in_si = self.units
        base = value.unit.base
        factor = value.unit.factor
        number = value.number
        new_unit = 1
        
        if not base:
            return value
            
        for n, unit in base:
            unit_in_generic, unit_in_si = self.find_generic_unit_for(unit)
            if not unit_in_si is None:
                factor = factor / (unit_in_si.factor ** n)
                new_unit *= (unit_in_generic.base[0][1] ** n)
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

class noconvert_generic_to_si(object):

    def __init__(self):
        pass

    def to_si(self, value):
        return value

    def to_generic(self, value):
        return value
