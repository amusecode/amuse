import numpy
from amuse.units.quantities import new_quantity
from amuse.units.quantities import is_unit
from amuse.units.quantities import is_quantity
from amuse.units.generic_unit_system import *

from amuse.support import exceptions
from amuse.support.core import late

class UnitsNotOrtogonalException(exceptions.AmuseException):
    formatstring = 'The number of orthoganal units is incorrect, expected {0} but found {1}. To convert between S.I. units and another system of units a set of quantities with orthogonal units is needed. These can be quantities with a single unit (such as length or time) or quantities with a derived units (such as velocity or force)'

class NotAQuantityException(exceptions.AmuseException):
    formatstring = 'Converters need to be initialized with a quantity argument[{0}] {1!r} is not a quantity'

class NotAScalarException(exceptions.AmuseException):
    formatstring = 'Converters need to be initialized with scalar quantities, argument[{0}] {1!r} is not a scalar'
    
class GenericToSiConverter(object):
    def __init__(self, generic_to_si):
        self.generic_to_si = generic_to_si
    
    def from_source_to_target(self, quantity):
        if hasattr(quantity, 'unit') and not quantity.unit.is_non_numeric():
            return self.generic_to_si.to_si(quantity) 
        else:
            return quantity
        
    def from_target_to_source(self, quantity):
        if hasattr(quantity, 'unit') and not quantity.unit.is_non_numeric():
            return self.generic_to_si.to_generic(quantity)
        else:
            return quantity
            
class SiToGenericConverter(object):
    def __init__(self, generic_to_si):
        self.generic_to_si = generic_to_si
    
    def from_source_to_target(self, quantity):
        if hasattr(quantity, 'unit') and not quantity.unit.is_non_numeric():
            return self.generic_to_si.to_generic(quantity) 
        else:
            return quantity
        
    def from_target_to_source(self, quantity):
        if hasattr(quantity, 'unit') and not quantity.unit.is_non_numeric():
            return self.generic_to_si.to_si(quantity)
        else:
            return quantity
                

class ConvertBetweenGenericAndSiUnits(object):
    """
    A ConvertBetweenGenericAndSiUnits object is a converter from
    arbitrary units which user gives, to si units (and vice versa).

    The ``generic_unit_converter'' **ConvertBetweenGenericAndSiUnits**
    is the actual class through which you define
    the unit system. Upon instantiation you choose the base units. In the
    example below we chose the speed of light as a unit: *c* = 1 unit length/second,
    and the second as the unit of time.

    Note that the system has two base dimensions, length and time. By the second argument we have
    assigned the unit second to time and by the requirement that unit lenght / second equals one,
    the new unit length will be {*c*} meters in S.I. units.

    Example::

        >>> from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
        >>> from amuse.units import units, constants
        >>> converter = ConvertBetweenGenericAndSiUnits(constants.c, units.m)

    """

    def __init__(self, *arguments_list):
        
        self.check_arguments(arguments_list)
        
        self.values = arguments_list
        
        self.units_of_values = [x.unit for x in self.values]
        self.system_rank = len(self.values)
        
        self.new_base = numpy.mat(numpy.zeros((self.system_rank,self.system_rank)))
        self.new_base_inv = numpy.mat(numpy.zeros((self.system_rank,self.system_rank)))
    
        available_units = set()
        for unit in self.units_of_values:
            for factor, base_unit in unit.base:
                available_units.add(base_unit)
        if not len(available_units) is self.system_rank:
            raise UnitsNotOrtogonalException(self.system_rank, len(available_units))
        self.list_of_available_units = list(available_units)
    
        self.new_base = self.determine_new_base()
        rank_of_new_base = self.matrixrank(self.new_base)
        if rank_of_new_base < self.system_rank:
            raise UnitsNotOrtogonalException(self.system_rank, rank_of_new_base)
        self.new_base_inv = self.new_base ** -1

    def check_arguments(self, arguments):
        
        for index, x in enumerate(arguments):
            if is_unit(x):
                continue
            if not is_quantity(x):
                raise NotAQuantityException(index, x)
            if not x.is_scalar():
                raise NotAScalarException(index, x)
        
    def matrixrank(self, A, tol=1e-8):
        s = numpy.linalg.svd(A,compute_uv=0)
        return numpy.sum(numpy.where( s>tol, 1, 0 ) )

    def determine_new_base(self):
        matrix = numpy.asmatrix(numpy.zeros((self.system_rank,self.system_rank)))
        for row, value in enumerate(self.values):
            for n, unit in value.unit.base:
                matrix[row, [i for i, j in enumerate(self.list_of_available_units) if j == unit]] = n
        return matrix

    @late
    def conversion_factors(self):
        factors_of_the_bases = numpy.asmatrix(numpy.zeros((self.system_rank,1)))
        for row, value in enumerate(self.values):
            factors_of_the_bases[row] = value.number * value.unit.factor
        log_factors_of_the_bases = numpy.log(numpy.abs(factors_of_the_bases))
        result = numpy.array(numpy.exp(self.new_base_inv*log_factors_of_the_bases))[:,0]
        return result

    @late
    def units(self):
        conversion_factors = self.conversion_factors
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
        for unit_generic, unit_in_si, base_in_si in self.unit_gerenic_to_unit_in_si:
            if base_in_si == unit:
                return unit_generic, unit_in_si
        return None, None


    def to_si(self, value):
        """     
        Convert a quantity in generic units to a quantity in
        S.I. units.
        
        .. code-block:: python
                          
            >>> from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
            >>> from amuse.units import units, constants
            >>> converter = ConvertBetweenGenericAndSiUnits(constants.c, units.s)
            >>> print converter.to_si(length)
            299792458.0 m
        """ 
        if not hasattr(value, 'unit'):
            return value
            
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

    def to_nbody(self, value):
        return self.to_generic(value)
        
    def to_generic(self, value):
        """
        Convert a quantity in S.I units to a quantity in
        generic units.
        
        .. code-block:: python
        
            >>> from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
            >>> from amuse.units import units, constants
            >>> converter = ConvertBetweenGenericAndSiUnits(constants.c, units.s)
            >>> print converter.to_generic(constants.c)
            1.0 length * time**-1

        """   
    
        generic_units_in_si = self.units
        if value.unit is None:
            return value
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

    def as_converter_from_si_to_generic(self):
        return SiToGenericConverter(self)

    def as_converter_from_generic_to_si(self):
        return GenericToSiConverter(self)

#    def __repr__(self):                                                                
#        string = ""                                                                    
#        string+= "generic unit converter object information"                           
#        string+= "-----------------------------------------"                           
#        string+= "unit system defined by setting {0} to ONE".format(self.values)       
#        string+= "Q = {Q}*[Q]; Quantity = value * unit"                                
#                                                                                       
#        for generic_unit, si_unit in self.units:                                       
#            string+= "Q = {{1}} * [{0}] Q = {{ {1} }} * [{2}^{3}]".format(             
#                generic_unit.base,                                                     
#                si_unit.number,                                                        
#                si_unit.base[1],                                                       
#                si_unit.base[0])                                                       
#                                                                                       
#        string+= self.list_of_available_units                                          
#        string+= self.new_base                                                         
#        string+= self.new_base_inv                                                     
#                                                                                       
#        return string                                                                  


    @late
    def mapping_from_base_in_si_to_unit_generic(self):
        result = {}
        for unit_generic, unit_in_si in self.units:
            base_in_si = unit_in_si.base[0][1]
            result[base_in_si] = unit_generic, unit_in_si      
        return result


    @late
    def unit_gerenic_to_unit_in_si(self):
        result = []
        for unit_generic, unit_in_si in self.units:
            base_in_si = unit_in_si.base[0][1]
            result.append((unit_generic, unit_in_si, base_in_si,))          
        return result


