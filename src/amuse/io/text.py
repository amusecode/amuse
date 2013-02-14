from amuse.support.core import late

from amuse.io import base
from amuse.units import units
from amuse.units import core
from amuse.units.quantities import is_quantity

import re

from amuse import datamodel
class LineBasedFileCursor(object):
    
    def __init__(self, file):
        self._line = None
        self.file = file
        self.read_next_line()
        
    def read_next_line(self):
        try:
            line = self.file.next()
            self._line = line.rstrip('\r\n')
        except StopIteration:
            self._line = None

    def forward(self):
        if not self.is_at_end():
            self.read_next_line()
    
    def line(self):
        return self._line
        
    def is_at_end(self):
        return self._line is None
        
        
class TableFormattedText(base.FileFormatProcessor):
    """
    Process a text file containing a table of values separated by a predefined character
    """
    
    provided_formats = ['txt']
    
    def __init__(self, filename = None, stream = None, set = None, format = None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        
        self.filename = filename
        self.stream = stream
        self.set = set
    
    def forward(self):
        line = self.data_file.readline()
        return line.rstrip('\r\n')
        
    def split_into_columns(self, line):
        if self.column_separator == ' ':
            return line.split()
        else:
            return line.split(self.column_separator)
    
    def load(self):
        if self.stream is None:
            self.stream = open(self.filename, "r")
            close_function = self.stream.close  
        else:
            close_function = lambda : None
            
        try:
            return self.load_from_stream()
        finally:
            close_function()
        
    def load_from_stream(self):
        self.cursor = LineBasedFileCursor(self.stream)
        self.read_header()
        self.read_rows()
        self.read_footer()
        return self.set
        
    def store(self):
        self.open_stream()
        try:
            return self.store_on_stream()
        finally:
            self.close_stream()
    
    def open_stream(self):
        if self.stream is None:
            self.stream = open(self.filename, "w")
            self.close_function = self.stream.close 
        else:
            self.close_function = lambda : None
            
    def close_stream(self):
        self.close_function()
        
    def store_on_stream(self):
        self.write_header()
        self.write_rows()
        self.write_footer()
        
    @base.format_option
    def attribute_names(self):
        """
        List of the names of the attributes to load or store
        """
        return self._get_attribute_names()
    
    def _get_attribute_names(self):
        if self.set is None:
            return map(lambda x : "col({0})".format(x), range(len(self.quantities)))
        else:
            return sorted(self.set.get_attribute_names_defined_in_store())
        
    @base.format_option
    def attribute_types(self):
        """
        List of the units of the attributes to store. If not given the 
        units will be derived from the units stored in the attribute set.
        When derived from the set, the units will always be converted to
        the base units.
        """
        return self._get_attribute_types()
    
    def _get_attribute_types(self):
        quantities = self.quantities
        if self.quantities:
            return map(lambda x : x.unit.to_simple_form() if is_quantity(x) else None, quantities)
        elif self.set is None:
            return [None] * len(self.attribute_names)
    
    @base.format_option
    def header_prefix_string(self):
        """
        Lines starting with this character will be handled as part of the header
        """
        return '#'
        
    
    @base.format_option
    def column_separator(self):
        """
        Separator between the columns
        """
        return ' '
        
    @base.format_option
    def footer_prefix_string(self):
        """
        Lines starting with this character will be handled as part of the footer
        """
        return self.header_prefix_string
        
    @base.format_option
    def key_in_column(self):
        """Column for the key (default is -1, no key stored/loaded).
        Keys will be interleaved with the data if set other than 0.
        No attribute type or name can be given for the key.
        """
        return -1
        
    
    @base.format_option
    def must_store_units_in_header(self):
        """Store line with parsable units in the header of the text file. This line can be used
        to restore the units on reading the file. When this line is added to the file you do not
        need to specify the attribute_types on reading the file"""
        return True
        
    def read_header(self):
        while not self.cursor.is_at_end() and self.cursor.line().startswith(self.header_prefix_string):
            self.read_header_line(self.cursor.line()[len(self.header_prefix_string):])
            self.cursor.forward()
    
    def read_header_line(self, line):
        pass
        
    def read_rows(self):
        values = map(lambda x : [], range(len(self.attribute_names)))
        
        number_of_particles = 0
        keys = []
        
        while not self.cursor.is_at_end() and not self.cursor.line().startswith(self.footer_prefix_string):
            columns = self.split_into_columns(self.cursor.line())
            if len(columns)>0:
                
                if self.key_in_column >= 0:
                       
                    if len(columns) != len(self.attribute_names) + 1:
                        raise base.IoException(
                            "Number of values on line '{0}' is {1}, expected {2}".format(self.cursor.line(), len(columns), len(self.attribute_names)))
            
                    
                    key = self.convert_string_to_long(columns[self.key_in_column])
                    keys.append(key)
                    del columns[self.key_in_column]
                    
                if len(columns) != len(self.attribute_names):
                    raise base.IoException(
                        "Number of values on line '{0}' is {1}, expected {2}".format(self.cursor.line(), len(columns), len(self.attribute_names)))
            
                for value_string, list_of_values in zip(columns, values):
                    list_of_values.append(self.convert_string_to_number(value_string))
                
                
                number_of_particles += 1
            self.cursor.forward()
    
        quantities = map(
            lambda value, unit : unit.new_quantity(value) if not unit is None else value, 
            values, 
            self.attribute_types
        )
        self.set = self.new_set(number_of_particles, keys = keys)
        self.set.set_values_in_store(self.set.get_all_indices_in_store(), self.attribute_names, quantities)
        
        self.cursor.forward()
        
    def read_footer(self):
        while not self.cursor.is_at_end() and self.cursor.line().startswith(self.footer_prefix_string):
            self.read_footer_line(self.cursor.line()[len(self.footer_prefix_string):])
            self.cursor.forward()
    
    def read_footer_line(self, line):
        pass
        
    def write_header(self):
        for x in self.header_lines():
            self.stream.write(self.header_prefix_string)
            self.stream.write(x)
            self.stream.write('\n')
        
    def write_rows(self):
        quantities = self.quantities
        units = self.attribute_types
        numbers = map(lambda quantity, unit : quantity if unit is None else quantity.value_in(unit), quantities, units)
        
        columns = []
        
        for x in numbers:
            columns.append(map(self.convert_number_to_string, x))
        
        rows = []
        for i in range(len(columns[0])):
            row = [x[i] for x in columns]
            
            if self.key_in_column >= 0:
                row.insert(self.key_in_column, self.convert_long_to_string(self.keys[i]))
            
            rows.append(row)
        
        lines = map(lambda  x : self.column_separator.join(x), rows)
        
        for x in lines:
            self.stream.write(x)
            self.stream.write('\n')
    
    def write_row(self, row):
        units = self.attribute_types
        row = map(lambda quantity, unit : quantity if unit is None else quantity.value_in(unit), row, units)
        row = map(self.convert_number_to_string, row)
        line = self.column_separator.join(row)
        self.stream.write(line)
        self.stream.write('\n')
        self.stream.flush()
        
        
    def write_footer(self):
        for x in self.footer_lines():
            self.stream.write(self.footer_prefix_string)
            self.stream.write(x)
            self.stream.write('\n')
        
    def header_lines(self):
        result = []
        if len(self.attribute_names) > 0:
            result.append(self.column_separator.join(self.attribute_names))
        result.append(self.column_separator.join(self.units_row))
        return result
        
    def footer_lines(self):
        return []
        
    def convert_string_to_number(self, string):
        return float(string)
        
    def convert_string_to_long(self, string):
        return long(string)
        
    def new_set(self, number_of_items, keys = []):
        if len(keys) > 0:
            return datamodel.Particles(number_of_items, keys = keys)
        else:
            return datamodel.Particles(number_of_items)
        
    @late
    def units_row(self):
        result = []
        for x in self.attribute_types:
            if x is None:
                result.append('-')
            else:
                result.append(str(x))
        return result
        
    @base.format_option
    def quantities(self):
        """
        List of vector quantities, each vector quantity is one column in the text file.
        By default this list will be derived from the particles set. When this
        option is given it will override the particle set data.
        """
        if self.set is None:
            return []
        else:
            return map(lambda x:getattr(self.set, x),self.attribute_names)


    @late
    def keys(self):
        if self.set is None:
            return []
        else:
            return self.set.key



        

    def convert_number_to_string(self, number):
        return str(number)
    

    def convert_long_to_string(self, number):
        return str(number)
    
    

    

    @base.format_option
    def float_format_string(self):
        "format specification string to convert numbers to strings, see format_spec in python documentation"
        return ".{0}e".format(self.precision_of_number_output)
    
    

    @base.format_option
    def precision_of_number_output(self):
        "The precision is a decimal number indicating how many digits should be displayed after the decimal point"
        return 12
    
    
class CsvFileText(TableFormattedText):
    """Process comma separated files
    
    Can process test files with comma separated fields.
    """
    
    provided_formats = ['csv']
    
    def __init__(self, filename = None, stream = None, set = None, format = None):
        TableFormattedText.__init__(self, filename, stream, set, format)
        if self.set is None:
            self.attribute_names = None
            self.attribute_types = None
    
    @base.format_option
    def column_separator(self):
        "separator between the columns"
        return ','
    

    def convert_string_to_unit(self, unit_string):
        if unit_string == '-':
            return None
        else:
            return eval(unit_string, core.__dict__)

    def convert_csv_string_to_unit(self, csv_string):
        return [self.convert_string_to_unit(sub) for sub in csv_string.split(self.column_separator)]
        
    def read_header_line(self, line):
        if self.attribute_names:
            if self.attribute_types is None:
                self.attribute_types = self.convert_csv_string_to_unit(line)
        else:
            self.attribute_names = [sub.strip() for sub in line.split(self.column_separator)]
    
    def header_lines(self):
        if self.attribute_names is None:
            self.attribute_names = self._get_attribute_names()
        if self.attribute_types is None:
            self.attribute_types = self._get_attribute_types()
        result = []
        result.append(self.column_separator.join(self.attribute_names))
        if self.must_store_units_in_header:
            result.append(self.column_separator.join(['-' if one_unit is None else one_unit.to_simple_form().reference_string() for one_unit in self.attribute_types]))
        result.append(self.column_separator.join(map(str, self.attribute_types)))
        return result
    
class Athena3DText(TableFormattedText):
    
    
    def read_header_line(self, line):
        line = line.lstrip()
        if line.startswith('['):
            self.read_attribute_names_from_line(line)
            
    def read_attribute_names_from_line(self, line):
        column_definition_strings = line.split()
        mapping_from_column_index_to_name = {}
        definition_re = re.compile(r'\[(\d+)\]=(.+)')
        for x in column_definition_strings:
            match_object = definition_re.match(x)
            if match_object is None:
                return
            index, name = match_object.groups()
            index = int(index) - 1
            name = name.replace('-','_')
            mapping_from_column_index_to_name[index] = name
        
        self.attribute_names = self.convert_dictionary_to_array(mapping_from_column_index_to_name)
        
    def convert_dictionary_to_array(self, dictionary):
        result = [None] * len(dictionary)
        for key, value in dictionary.iteritems():
            result[key] = value
        
        return result
        
            
        
    
