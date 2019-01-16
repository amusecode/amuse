import numpy
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
        if not stream is None:
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
    def attribute_dtypes(self):
        """
        List of the data types of the attributes. If not given, the 
        data types will be float64.
        """
        return self._get_attribute_dtypes()
    
    def _get_attribute_dtypes(self):
        return [None] * len(self.attribute_names)
    
    def _new_converter_from_string_to_dtype(self, dtype):
        if dtype is None:
            return self.convert_string_to_number
        kind = numpy.array([], dtype=dtype).dtype.kind
        if kind == 'f':
            return self.convert_string_to_number
        elif kind == 'S' or kind == 'U':
            return lambda string_value: string_value
        elif kind == 'b':
            return lambda string_value: string_value == 'True' or string_value == 'true'
        else:
            return self.convert_string_to_long
    

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
            header_line = self.cursor.line()[len(self.header_prefix_string):]
            if not header_line.startswith(self.header_prefix_string):
                self.read_header_line(header_line)
            self.cursor.forward()
    

    def read_header_line(self, line):
        line = line.strip()
        if line == 'AMUSE CSV: 1.0':
            self.has_amuse_header = True
        if self.has_amuse_header:
            KEY_IN_COLUMN_HDR='KEY IN COLUMN:'
            if line.startswith(KEY_IN_COLUMN_HDR):
                self.key_in_column = int(line[len(KEY_IN_COLUMN_HDR):].strip())
            
            COLUMN_HDR='COL:'
            if line.startswith(COLUMN_HDR):
                line = line[len(COLUMN_HDR):].strip()
                parts = line.split(':')
                parts = map(lambda x : x.strip(), parts)
                index = int(parts[0])
                if self.key_in_column >= 0 and index > self.key_in_column:
                    index -= 1
                name = parts[1]
                is_on = parts[2] == 'on'
                dtype = map(float, parts[3:-1])
                unit = self.convert_float_to_unit(dtype)
                if self.attribute_names is None:
                    current = []
                else:
                    current = self.attribute_names[:]
                print index
                while len(current)<=index:
                    current.append("")
                current[index] = name
                self.attribute_names = current
                if self.attribute_types is None:
                    current = []
                else:
                    current = self.attribute_types[:]
                while len(current)<=index:
                    current.append(None)
                current[index] = unit
                self.attribute_types = current
                self.attribute_dtypes = self._get_attribute_dtypes()
        

    def read_rows(self):
        values = map(lambda x : [], range(len(self.attribute_names)))
        
        number_of_particles = 0
        keys = []
        
        self.set = None
        string_converters = map(self._new_converter_from_string_to_dtype, self.attribute_dtypes)
        units_with_dtype = map(core.unit_with_specific_dtype, self.attribute_types, self.attribute_dtypes)
        
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
            
                map(lambda value_string, list_of_values, conv: list_of_values.append(conv(value_string)), 
                    columns, values, string_converters)
                
                number_of_particles += 1
            self.cursor.forward()
            if number_of_particles >= self.maximum_number_of_lines_buffered:
                quantities = map(
                    lambda value, unit, dtype : unit.new_quantity(value) if not unit is None else (numpy.asarray(value, dtype=dtype) if not dtype is None else value), 
                    values, 
                    units_with_dtype,
                    self.attribute_dtypes
                )
                if self.set is None:
                    self.set = self.new_set(number_of_particles, keys = keys)
                    self.set.set_values_in_store(self.set.get_all_indices_in_store(), self.attribute_names, quantities)
                else:
                    tmp_set = self.new_set(number_of_particles, keys = keys)
                    tmp_set.set_values_in_store(tmp_set.get_all_indices_in_store(), self.attribute_names, quantities)
                    self.set.add_particles(tmp_set)
                    
                number_of_particles = 0
                keys = []
                values = map(lambda x : [], range(len(self.attribute_names)))
        
        if number_of_particles > 0:
            quantities = map(
                lambda value, unit, dtype : unit.new_quantity(value) if not unit is None else (numpy.asarray(value, dtype=dtype) if not dtype is None else value), 
                values, 
                units_with_dtype,
                self.attribute_dtypes
            )
            if self.set is None:
                self.set = self.new_set(number_of_particles, keys = keys)
                self.set.set_values_in_store(self.set.get_all_indices_in_store(), self.attribute_names, quantities)
            else:
                tmp_set = self.new_set(number_of_particles, keys = keys)
                tmp_set.set_values_in_store(tmp_set.get_all_indices_in_store(), self.attribute_names, quantities)
                self.set.add_particles(tmp_set)
        elif self.set is None:
            self.set = self.new_set(0)
            
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
        
        max_row = -1
        for x in quantities:
            max_row = max(max_row, len(x))
        
        block_size = min(self.maximum_number_of_lines_buffered, max_row)
        offset = 0
        
        while offset < max_row:
            
            numbers = map(lambda quantity, unit : quantity[offset:offset+block_size] if unit is None else quantity[offset:offset+block_size].value_in(unit), quantities, units)
            
            columns = []
            
            
                
            for x in numbers:
                columns.append(map(self.convert_number_to_string, x))
            
            rows = []
            for i in range(len(columns[0])):
                row = [x[i] for x in columns]
                
                if self.key_in_column >= 0:
                    row.insert(self.key_in_column, self.convert_long_to_string(self.keys[i+offset]))
                
                rows.append(row)
            
            lines = map(lambda  x : self.column_separator.join(x), rows)
            
            
            for x in lines:
                self.stream.write(x)
                self.stream.write('\n')
                
            offset += block_size
        
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
        if self.has_amuse_header:
            result.append('AMUSE CSV: 1.0')
            result.append('KEY IN COLUMN: '+str(self.key_in_column))
            result.append('HEADERS:')
            for i, (x, unit) in enumerate(zip(self.attribute_names, self.attribute_types)):
                column = i
                if self.key_in_column >= 0:
                    if i >= self.key_in_column:
                        column += 1
                if unit is None:
                    unitstring= '1:-1:0:0:0:0:0:0:0'
                    description = '(-)'
                else:
                    unitstring = '{0:.18g}:{1:.0f}:{2:.18g}:{3:.18g}:{4:.18g}:{5:.18g}:{6:.18g}:{7:.18g}:{8:.18g}'.format(*unit.to_array_of_floats())
                    description = '(' + unit.describe_array_of_floats() + ')'
                result.append('COL:{0}:{1}:{2}:{3}:{4}'.format(column, x, 'on', unitstring, description))
            
            result.append('')
        else:
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
            if len(self.set) == 0:
                return map(lambda x:[],self.attribute_names)
            else:
                return map(lambda x:getattr(self.set, x),self.attribute_names)


    @late
    def keys(self):
        if self.set is None:
            return []
        else:
            return self.set.key



        

    def convert_number_to_string(self, number):
        if self.is_precise:
            return '{:.18e}'.format(number)
        else:
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
    
    @base.format_option
    def maximum_number_of_lines_buffered(self):
        """"The maximum number of lines to convert and write out in one go, larger values may be a little faster but will
        use a lot more memory"""
        return 10000
    
    
    
    @base.format_option
    def is_precise(self):
        """"Output most precise text output for floats (9 digitst) and doubles (17 digits)"""
        return False




    @base.format_option
    def has_amuse_header(self):
        """If true, store the decription of the attributes and units in the header, better format than must_store_units_in_header provides """
        return False


    @base.format_option
    def comments(self):
        """
        Comments to add to the header
        """
        return []




    def convert_float_to_unit(self, floats):
        from amuse.units import core
        from amuse.units import units
        
        print floats, int(floats[1]) == -1, numpy.all(numpy.asarray(floats[2:]) == 0.0), numpy.asarray(floats[2:]) == 0.0
        if int(floats[1]) == -1 and numpy.all(numpy.asarray(floats[2:]) == 0.0):
            return None
        factor = floats[0]
        result = factor
        system_index = int(floats[1])
        unit_system = None
        for x in core.system.ALL.values():
            if x.index == system_index:
                unit_system = x
                break
        for x in unit_system.bases:
            power = floats[x.index + 2]
            if not power == 0.0:
                result = result * (x ** power)
        return result
        








    @base.format_option
    def stream(self):
        """"Set the stream to output to"""
        return None





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
            try:
                return eval(unit_string, core.__dict__)
            except:
                if hasattr(units, unit_string):
                    return getattr(units, unit_string)
                return eval(unit_string, units.__dict__)

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
        if len(self.comments) > 0:
            result.extend(map(lambda x : self.header_prefix_string + x, self.comments))
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
        
            
        
    
class AmuseText(TableFormattedText):
    """Process text files with AMUSE header format
    """
    
    provided_formats = ['amuse-txt']
    
    
    @base.format_option
    def has_amuse_header(self):
        """If true, store the decription of the attributes and units in the header, better format than must_store_units_in_header provides (default True for amuse csv format)"""
        return True







    @base.format_option
    def is_precise(self):
        """"Output most precise text output for floats (9 digitst) and doubles (17 digits) (True for amuse text)"""
        return True





    def convert_number_to_string(self, number):
        if self.is_precise:
            return '{:.18e}'.format(number)
        else:
            return str(number)

    @base.format_option
    def use_fractions(self):
        """"Output floats as fractions, will be more precise"""
        return True
