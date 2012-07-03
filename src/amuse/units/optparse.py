from __future__ import absolute_import

import optparse
import textwrap
from amuse.units import quantities

def check_builtin_unit(option, opt, value):
    (cvt, what) = optparse._builtin_cvt[option.type]
    try:
        result = cvt(value)
        if option.unit is None:
            return result
        else:
            return quantities.new_quantity(result, option.unit)
            
    except ValueError:
        raise OptionValueError(
            _("option %s: invalid %s value: %r") % (opt, what, value))


class Option(optparse.Option):
    TYPE_CHECKER = {
        "int"    : check_builtin_unit,
        "long"   : check_builtin_unit,
        "float"  : check_builtin_unit,
        "complex": check_builtin_unit,
        "choice" : optparse.check_choice,
    }
    ATTRS = optparse.Option.ATTRS + ['unit',]
    
    def convert_value(self, opt, value):
        if value is not None:
            if self.nargs == 1:
                return self.check_value(opt, value)
            else:
                return tuple([self.check_value(opt, v) for v in value])

class IndentedHelpFormatter(optparse.IndentedHelpFormatter):
    
    def __init__(self,
             indent_increment=2,
             max_help_position=24,
             width=None,
             short_first=1):
        optparse.IndentedHelpFormatter.__init__(
            self, indent_increment, max_help_position, width, short_first)
            
        self.unit_tag = "%unit"
        
    def format_option(self, option):
        # The help for each option consists of two parts:
        #   * the opt strings and metavars
        #     eg. ("-x", or "-fFILENAME, --file=FILENAME")
        #   * the user-supplied help string
        #     eg. ("turn on expert mode", "read data from FILENAME")
        #
        # If possible, we write both of these on the same line:
        #   -x      turn on expert mode
        #
        # But if the opt string list is too long, we put the help
        # string on a second line, indented to the same column it would
        # start in if it fit on the first line.
        #   -fFILENAME, --file=FILENAME
        #           read data from FILENAME
        result = []
        opts = self.option_strings[option]
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else:                       # start help on same line as opts
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0
        result.append(opts)
            
        if option.help:
            help_text = self.expand_default(option)
            help_text = self.expand_unit(option, help_text)
            help_lines = textwrap.wrap(help_text, self.help_width)
            result.append("%*s%s\n" % (indent_first, "", help_lines[0]))
            result.extend(["%*s%s\n" % (self.help_position, "", line)
                           for line in help_lines[1:]])
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)
    
    def expand_unit(self, option, help_text):
        if self.parser is None or not self.unit_tag:
            return help_text

        if option.unit is None:
            return help_text
        
        unit = option.unit

        return help_text.replace(self.unit_tag, str(unit))
    
    
    
class OptionParser(optparse.OptionParser):
    
    def __init__(self,
             usage=None,
             option_list=None,
             option_class=Option,
             version=None,
             conflict_handler="error",
             description=None,
             formatter=None,
             add_help_option=True,
             prog=None,
             epilog=None):
        
        
        if formatter is None:
            formatter = IndentedHelpFormatter()
            
        optparse.OptionParser.__init__(
            self, usage, 
            option_list, option_class, 
            version, conflict_handler, 
            description, formatter,
            add_help_option, prog, epilog
        )

    def get_default_values(self):
        if not self.process_default_values:
                # Old, pre-Optik 1.5 behaviour.
            return Values(self.defaults)

        defaults = self.defaults.copy()
        for option in self._get_all_options():
            default = defaults.get(option.dest)
            if optparse.isbasestring(default):
                opt_str = option.get_opt_string()
                defaults[option.dest] = option.check_value(opt_str, default)
            elif not option.unit is None and not quantities.is_quantity(default):
                defaults[option.dest] = quantities.new_quantity(default, option.unit)
                
        return optparse.Values(defaults)
