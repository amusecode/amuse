import numpy
from amuse.support.exceptions import AmuseException
from amuse.support import options

registered_printing_strategies = {}

class UnsupportedPrintingStrategyException(AmuseException):
    """
    Raised when the given printing strategy is not supported by AMUSE.
    """
    formatstring = "You tried to set the printing strategy to '{0}', but this printing strategy is not available"



class PrintingStrategy(object):
    
    def __init__(self):
        pass
    
    def convert_quantity(self, quantity):
        return quantity
    
    def quantity_to_string(self, quantity):
        quantity = self.convert_quantity(quantity)
        result = self.string_prefix()
        result += self.string_number(quantity)
        string_of_unit = self.string_unit(quantity)
        if string_of_unit:
            result += self.string_separator()
        result += string_of_unit
        result += self.string_suffix()
        return result
    
    def string_prefix(self):
        return ""
    
    def string_number(self, quantity):
        return ""
        
    def string_separator(self):
        return " "
    
    def string_unit(self, quantity):
        return ""
    
    def string_suffix(self):
        return ""
    
    def old_numbers_to_string(self, number):
        if hasattr(number, "__iter__"):
            return '[' + ', '.join([str(x) for x in number]) + ']'
        return str(number)
    
    def numbers_to_string(self, quantity, precision=None):
        if quantity.is_vector():
            if precision is None:
                return '[' + ', '.join(["%s" % x for x in quantity.number]) + ']'
            else:
                fmt = "%#." + str(precision) + "g"
                return '[' + ', '.join([fmt % x for x in quantity.number]) + ']'
        else:
            if precision is None:
                return "%s" % quantity.number
            else:
                return ("%."+str(precision)+"g") % quantity.number
    
    @classmethod
    def register(cls):
        """
        Register this class, so that it can be found by name
        in the :func:`set_printing_strategy` function.
        """
        add_printing_strategy(cls)
    


class DefaultPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['default', 'with_units']
    
    def string_number(self, quantity):
        return self.numbers_to_string(quantity)
    
    def string_unit(self, quantity):
        return str(quantity.unit)
    



class SimplePrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['simple']
    
    def string_number(self, quantity):
        factor, print_unit = quantity.unit.to_factor_and_reduced_form()
        return self.numbers_to_string(quantity * factor)
    
    def string_unit(self, quantity):
        factor, print_unit = quantity.unit.to_factor_and_reduced_form()
        return str(print_unit)
        
class NoUnitsPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['no_unit', 'no_units']
    
    def string_number(self, quantity):
        return self.numbers_to_string(quantity)
    


class FormalPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['formal',]
    
    def string_prefix(self):
        return "<quantity "
    
    def string_number(self, quantity):
        return self.numbers_to_string(quantity)
    
    def string_separator(self):
        return " | "
    
    def string_unit(self, quantity):
        return str(quantity.unit)
    
    def string_suffix(self):
        return ">"
    


class NBodyPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['nbody',]
    
    def __init__(self, nbody_converter = None, ignore_converter_exceptions = False):
        self.ignore_converter_exceptions = ignore_converter_exceptions
        self.nbody_converter = nbody_converter
    
    def convert_quantity(self, quantity):
        if is_not_nbody_unit(quantity.unit):
            if self.nbody_converter:
                return self.nbody_converter.to_nbody(quantity)
            elif not self.ignore_converter_exceptions:
                raise AmuseException("Unable to convert {0} to N-body units. No "
                    "nbody_converter given".format(quantity.unit))
        return quantity
    
    def string_number(self, quantity):
        return self.numbers_to_string(quantity)
    


class PrintingStrategyWithPreferredUnits(PrintingStrategy):
    
    def convert_quantity(self, quantity):
        if has_nbody_unit(quantity.unit):
            if self.nbody_converter:
                return _quantity_in_preferred_units(self.preferred_units, self.nbody_converter.to_si(quantity))
            elif not self.ignore_converter_exceptions:
                raise AmuseException("Unable to convert {0} to SI units. No "
                    "nbody_converter given".format(quantity.unit))
        return _quantity_in_preferred_units(self.preferred_units, quantity)
    
    def string_number(self, quantity):
        return self.numbers_to_string(quantity)
    
    def string_unit(self, quantity):
        if self.print_units:
            return str(quantity.unit)
        else:
            return ""
    


class AstroPrintingStrategy(PrintingStrategyWithPreferredUnits):
    
    provided_strategy_names = ['astro',]
    
    def __init__(self, nbody_converter = None, print_units = True, ignore_converter_exceptions = None):
        self.ignore_converter_exceptions = (print_units if (ignore_converter_exceptions is None) 
            else ignore_converter_exceptions)
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        from amuse.units import units
        self.preferred_units = [units.MSun, units.Myr, units.parsec, units.J]
    


class SIPrintingStrategy(PrintingStrategyWithPreferredUnits):
    
    provided_strategy_names = ['SI', 'si', 'MKS', 'mks']
    
    def __init__(self, nbody_converter = None, print_units = True, ignore_converter_exceptions = None):
        self.ignore_converter_exceptions = (print_units if (ignore_converter_exceptions is None) 
            else ignore_converter_exceptions)
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        from amuse.units import units
        self.preferred_units = [units.m, units.kg, units.s, units.A, units.K, units.mol, units.cd]
    


class CGSPrintingStrategy(PrintingStrategyWithPreferredUnits):
    
    provided_strategy_names = ['CGS', 'cgs']
    
    def __init__(self, nbody_converter = None, print_units = True, ignore_converter_exceptions = None):
        self.ignore_converter_exceptions = (print_units if (ignore_converter_exceptions is None) 
            else ignore_converter_exceptions)
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        from amuse.units import units
        self.preferred_units = [units.cm, units.g, units.s, units.A, units.K, units.mol, units.cd]
    


class CustomPrintingStrategy(PrintingStrategyWithPreferredUnits):
    
    provided_strategy_names = ['custom',]
    
    def __init__(self, nbody_converter=None, print_units=True, preferred_units=None,
            precision=None, prefix="", separator=" ", suffix="",
            ignore_converter_exceptions=None):
        self.ignore_converter_exceptions = (print_units if (ignore_converter_exceptions is None) 
            else ignore_converter_exceptions)
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        self.preferred_units = preferred_units
        self.precision = precision
        self.prefix = prefix
        self.separator = separator
        self.suffix = suffix
    
    def string_prefix(self):
        return self.prefix
    
    def string_number(self, quantity):
        return self.numbers_to_string(quantity, precision = self.precision)
    
    def string_separator(self):
        return self.separator
    
    def string_suffix(self):
        return self.suffix


def _quantity_in_preferred_units(preferred_units, quantity):
    if preferred_units is None:
        return quantity
    for preferred_unit in preferred_units:
        if quantity.unit.has_same_base_as(preferred_unit):
            return quantity.as_quantity_in(preferred_unit)
    if "factor_unit" in str(quantity.unit.__class__):
        local = _quantity_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0))
        return local.unit.new_quantity(quantity.number * quantity.unit.local_factor * local.number)
    if "mul_unit" in str(quantity.unit.__class__):
        left  = _quantity_in_preferred_units(preferred_units, quantity.unit.left_hand(1.0))
        right = _quantity_in_preferred_units(preferred_units, quantity.unit.right_hand(1.0))
        return (left.unit * right.unit).new_quantity(quantity.number * left.number * right.number)
    if "div_unit" in str(quantity.unit.__class__):
        left  = _quantity_in_preferred_units(preferred_units, quantity.unit.left_hand(1.0))
        right = _quantity_in_preferred_units(preferred_units, quantity.unit.right_hand(1.0))
        return (left.unit / right.unit).new_quantity(quantity.number * left.number / right.number)
    if "pow_unit" in str(quantity.unit.__class__):
        local = _quantity_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0))**quantity.unit.power
        return local.unit.new_quantity(quantity.number * local.number)
    if "named_unit" in str(quantity.unit.__class__):
        local = _quantity_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0))
        return local.unit.new_quantity(quantity.number * local.number)
    return quantity

def has_nbody_unit(unit):
    for factor, x in unit.base:
        if x.is_generic():
            return True
    return False

def is_not_nbody_unit(unit):
    for factor, x in unit.base:
        if not x.is_generic():
            return True
    return False


def set_printing_strategy(strategy, **kwargs):
    global current_printing_strategy
    current_printing_strategy = _get_printing_strategy_factory(strategy)(**kwargs)

def get_current_printing_strategy():
    global current_printing_strategy
    return current_printing_strategy.__class__


def _get_printing_strategy_factory(strategy):
    if isinstance(strategy, basestring):
        if not strategy in registered_printing_strategies:
            raise UnsupportedPrintingStrategyException(strategy)
        return registered_printing_strategies[strategy]
    else:
        return strategy

def add_printing_strategy(class_of_the_printing_strategy):
    """
    Register the specified class, so that it can be used
    by the :func:`set_printing_strategy` function.
    
    Do not call this method directly, instead use :func:`PrintingStrategy.register`
    """
    for x in class_of_the_printing_strategy.provided_strategy_names:
        registered_printing_strategies[x] = class_of_the_printing_strategy


class _Defaults(options.OptionalAttributes):
    
    @options.option(sections=['output',])
    def printing_strategy(self):
        return 'default'

DefaultPrintingStrategy.register()
NoUnitsPrintingStrategy.register()
FormalPrintingStrategy.register()
NBodyPrintingStrategy.register()
AstroPrintingStrategy.register()
SIPrintingStrategy.register()
CGSPrintingStrategy.register()
CustomPrintingStrategy.register()
SimplePrintingStrategy.register()

set_printing_strategy(_Defaults().printing_strategy)
