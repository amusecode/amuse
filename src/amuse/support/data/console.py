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
    


def set_printing_stategy(strategy):
    global current_printing_stategy
    current_printing_stategy = _get_printing_strategy_factory(strategy)()


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
    by the :func:`set_printing_stategy` function.
    
    Do not call this method directly, instead use :func:`PrintingStrategy.register`
    """
    for x in class_of_the_printing_strategy.provided_strategy_names:
        registered_printing_strategies[x] = class_of_the_printing_strategy


DefaultPrintingStrategy.register()
NoUnitsPrintingStrategy.register()

set_printing_stategy('default')
