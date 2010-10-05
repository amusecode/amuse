from amuse.support.exceptions import AmuseException


registered_printing_strategies = {}

class UnsupportedPrintingStrategyException(AmuseException):
    """
    Raised when the given printing strategy is not supported by AMUSE.
    """
    formatstring = "You tried to set the printing strategy to '{0}', but this printing strategy is not available"



class PrintingStrategy(object):
    
    def __init__(self):
        pass
    
    def quantity_to_string(self, quantity):
        result = self.string_prefix()
        result += self.string_number(quantity)
        result += self.string_separator()
        result += self.string_unit(quantity)
        result += self.string_suffix()
        return result
    
    def string_prefix(self):
        return ""
    
    def string_number(self, quantity):
        return ""
        
    def string_separator(self):
        return ""
    
    def string_unit(self, quantity):
        return ""
    
    def string_suffix(self):
        return ""
    
    def default_string_converter_for_numbers(self, number):
        if hasattr(number, "__iter__"):
            return '[' + ', '.join([str(x) for x in number]) + ']'
        return str(number)
    
    @classmethod
    def register(cls):
        """
        Register this class, so that it can be found by name
        int the :func:`write_set_to_file` and :func:`read_set_from_file`
        functions.
        """
        add_printing_strategy(cls)
    


class DefaultPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['default', 'with_units']
    
    def string_number(self, quantity):
        return self.default_string_converter_for_numbers(quantity.number)
    
    def string_unit(self, quantity):
        result = str(quantity.unit)
        if result:
            return " " + result
    


class NoUnitsPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['no_unit', 'no_units']
    
    def string_number(self, quantity):
        return self.default_string_converter_for_numbers(quantity.number)
    


class FormalPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['formal',]
    
    def string_prefix(self):
        return "<quantity "
    
    def string_number(self, quantity):
        return self.default_string_converter_for_numbers(quantity.number)
    
    def string_separator(self):
        return " | "
    
    def string_unit(self, quantity):
        return str(quantity.unit)
    
    def string_suffix(self):
        return ">"
    


class NBodyPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['nbody',]
    
    def __init__(self, nbody_converter = None):
        self.nbody_converter = nbody_converter
    
    def is_not_nbody_unit(self, unit):
        for factor, x in unit.base:
            if not x.is_generic():
                return True
        return False
    
    def string_number(self, quantity):
        if self.is_not_nbody_unit(quantity.unit):
            if self.nbody_converter:
                quantity = self.nbody_converter.to_nbody(quantity)
            else:
                raise AmuseException("Unable to convert {0} to N-body units. No "
                    "nbody_converter given".format(quantity.unit))
        return self.default_string_converter_for_numbers(quantity.number)
    


class AstroPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['astro',]
    
    def __init__(self, nbody_converter = None, print_units = True):
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        from amuse.support.units import units
        self.supported_units = [units.MSun, units.Myr, units.parsec, units.J]
    
    def has_nbody_unit(self, unit):
        for factor, x in unit.base:
            if x.is_generic():
                return True
        return False
    
    def string_number(self, quantity):
        if self.has_nbody_unit(quantity.unit):
            if self.nbody_converter:
                quantity = self.nbody_converter.to_si(quantity)
            else:
                raise AmuseException("Unable to convert {0} to SI units. No "
                    "nbody_converter given".format(quantity.unit))
        return self.default_string_converter_for_numbers(
            _number_in_preferred_units(self.supported_units, quantity))
    
    def string_unit(self, quantity):
        if self.print_units:
            if self.has_nbody_unit(quantity.unit):
                quantity = self.nbody_converter.to_si(quantity)
            return " " + str(_unit_in_preferred_units(self.supported_units, quantity))
        else:
            return ""
    


class SIPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['SI', 'si']
    
    def __init__(self, nbody_converter = None, print_units = True):
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        from amuse.support.units import units
        self.supported_units = [units.m, units.kg, units.s, units.A, units.K, units.mol, units.cd]
    
    def has_nbody_unit(self, unit):
        for factor, x in unit.base:
            if x.is_generic():
                return True
        return False
    
    def string_number(self, quantity):
        if self.has_nbody_unit(quantity.unit):
            if self.nbody_converter:
                quantity = self.nbody_converter.to_si(quantity)
            else:
                raise AmuseException("Unable to convert {0} to SI units. No "
                    "nbody_converter given".format(quantity.unit))
        return self.default_string_converter_for_numbers(
            _number_in_preferred_units(self.supported_units, quantity))
    
    def string_unit(self, quantity):
        if self.print_units:
            if self.has_nbody_unit(quantity.unit):
                quantity = self.nbody_converter.to_si(quantity)
            return " " + str(_unit_in_preferred_units(self.supported_units, quantity))
        else:
            return ""
    


class CustomPrintingStrategy(PrintingStrategy):
    
    provided_strategy_names = ['custom',]
    
    def __init__(self, nbody_converter = None, print_units = True, preferred_units = [], digits = None,
            prefix = "", separator = " ", suffix = ""):
        self.nbody_converter = nbody_converter
        self.print_units = print_units
        self.preferred_units = preferred_units
        self.digits = digits
        self.prefix = prefix
        self.separator = separator
        self.suffix = suffix
    
    def has_nbody_unit(self, unit):
        for factor, x in unit.base:
            if x.is_generic():
                return True
        return False
    
    def custom_string_converter_for_numbers(self, number, digits = None):
        if digits:
            format = "%." + str(digits) + "g"
        else:
            format = "%s"
        if hasattr(number, "__iter__"):
            return '[' + ', '.join([format % x for x in number]) + ']'
        return format % number
    
    def string_prefix(self):
        return self.prefix
    
    def string_number(self, quantity):
        if self.has_nbody_unit(quantity.unit):
            if self.nbody_converter:
                quantity = self.nbody_converter.to_si(quantity)
            else:
                raise AmuseException("Unable to convert {0} to SI units. No "
                    "nbody_converter given".format(quantity.unit))
        return self.custom_string_converter_for_numbers(
            _number_in_preferred_units(self.preferred_units, quantity), digits = self.digits)
    
    def string_separator(self):
        return self.separator
    
    def string_unit(self, quantity):
        if self.print_units:
            if self.has_nbody_unit(quantity.unit):
                quantity = self.nbody_converter.to_si(quantity)
            return str(_unit_in_preferred_units(self.preferred_units, quantity))
        else:
            return ""
    
    def string_suffix(self):
        return self.suffix


def _number_in_preferred_units(preferred_units, quantity):
    if quantity.unit in preferred_units:
        return quantity.number
    for supported_unit in preferred_units:
        if quantity.unit.are_bases_equal(supported_unit):
            return quantity.as_quantity_in(supported_unit).number
    if "factor_unit" in str(quantity.unit.__class__):
        return (quantity.number * quantity.unit.local_factor * 
            _number_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0)))
    if "mul_unit" in str(quantity.unit.__class__):
        return (quantity.number * 
            _number_in_preferred_units(preferred_units, quantity.unit.left_hand(1.0)) * 
            _number_in_preferred_units(preferred_units, quantity.unit.right_hand(1.0)))
    if "div_unit" in str(quantity.unit.__class__):
        return (quantity.number * 
            _number_in_preferred_units(preferred_units, quantity.unit.left_hand(1.0)) / 
            _number_in_preferred_units(preferred_units, quantity.unit.right_hand(1.0)))
    if "pow_unit" in str(quantity.unit.__class__):
        return quantity.number * _number_in_preferred_units(preferred_units, 
            quantity.unit.local_unit(1.0))**quantity.unit.power
    if "named_unit" in str(quantity.unit.__class__):
        return quantity.number * _number_in_preferred_units(preferred_units, 
            quantity.unit.local_unit(1.0))
    return quantity.number

def _unit_in_preferred_units(preferred_units, quantity):
    if quantity.unit in preferred_units:
        return quantity.unit
    for supported_unit in preferred_units:
        if quantity.unit.are_bases_equal(supported_unit):
            return supported_unit
    if "factor_unit" in str(quantity.unit.__class__):
        return _unit_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0))
    if "mul_unit" in str(quantity.unit.__class__):
        return (_unit_in_preferred_units(preferred_units, quantity.unit.left_hand(1.0)) * 
            _unit_in_preferred_units(preferred_units, quantity.unit.right_hand(1.0)))
    if "div_unit" in str(quantity.unit.__class__):
        return (_unit_in_preferred_units(preferred_units, quantity.unit.left_hand(1.0)) / 
            _unit_in_preferred_units(preferred_units, quantity.unit.right_hand(1.0)))
    if "pow_unit" in str(quantity.unit.__class__):
        return _unit_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0))**quantity.unit.power
    if "named_unit" in str(quantity.unit.__class__):
        return _unit_in_preferred_units(preferred_units, quantity.unit.local_unit(1.0))
    return quantity.unit


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


DefaultPrintingStrategy.register()
NoUnitsPrintingStrategy.register()
FormalPrintingStrategy.register()
NBodyPrintingStrategy.register()
AstroPrintingStrategy.register()
SIPrintingStrategy.register()
CustomPrintingStrategy.register()

set_printing_strategy('default')
