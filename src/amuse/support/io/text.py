from amuse.support.data import core
from amuse.support.core import late
from amuse.support.units import units

import re

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
        
        
class TableFormattedText(object):
    
    def __init__(self, filename, stream = None, set = None):
        self.filename = filename
        self.stream = stream
        self.set = set
    
    def forward(self):
        line = selfdata_file.readline()
        return line.rstrip('\r\n')
        
    def split_into_columns(self, line):
        return line.split()
    
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
        
        if self.stream is None:
            self.stream = open(self.filename, "w")
            close_function = self.stream.close 
        else:
            close_function = lambda : None
            
        try:
            return self.store_on_stream()
        finally:
            close_function()
            
    def store_on_stream(self):
        self.write_header()
        self.write_rows()
        self.write_footer()
        
    @late
    def attribute_names(self):
        if self.set is None:
            return []
        else:
            return self.set.stored_attributes()
        
    @late
    def attribute_types(self):
        if self.set is None:
            return map(lambda x : units.none , self.attribute_names)
        else:
            return map(lambda x : units.none , self.attribute_names)
    
    @late
    def header_prefix_string(self):
        return '#'
        
    
    @late
    def column_separator(self):
        return ' '
        
    @late
    def footer_prefix_string(self):
        return self.header_prefix_string
        
    def read_header(self):
        while not self.cursor.is_at_end() and self.cursor.line().startswith(self.header_prefix_string):
            self.read_header_line(self.cursor.line()[len(self.header_prefix_string):])
            self.cursor.forward()
    
    def read_header_line(self, line):
        pass
        
    def read_rows(self):
        values = map(lambda x : [], range(len(self.attribute_names)))
        
        number_of_particles = 0
        
        while not self.cursor.is_at_end() and not self.cursor.line().startswith(self.footer_prefix_string):
            columns = self.split_into_columns(self.cursor.line())
            if len(columns) != len(self.attribute_names):
                raise Exception(
                    "Number of values on line '{0}' is {1}, expected {2}".format(self.cursor.line(), len(columns), len(self.attribute_names)))
            
            
            for value_string, list_of_values in zip(columns, values):
                list_of_values.append(self.convert_string_to_number(value_string))
                
            self.cursor.forward()
            number_of_particles += 1
        
        quantities = map(
            lambda value, unit : unit.new_quantity(value), 
            values, 
            self.attribute_types
        )
        self.set = self.new_set(number_of_particles)
        self.set._set_values(self.set._get_keys(), self.attribute_names, quantities)
        
        self.cursor.forward()
        
    def read_footer(self):
        while not self.cursor.is_at_end() and self.cursor.line().startswith(self.footer_prefix_string):
            self.read_footer_line(self.line()[len(self.footer_prefix_string):])
            self.cursor.forward()
    
    def read_footer_line(self, line):
        pass
        
    def write_header(self):
        for x in self.header_lines():
            self.stream.write(self.header_prefix_string)
            self.stream.write(x)
            self.stream.write('\n')
        
    def write_rows(self):
        keys = self.set.key
        quantities = map(lambda x:getattr(self.set, x),self.attribute_names)
        units = self.attribute_types
        numbers = map(lambda quantity, unit : quantity.value_in(unit), quantities, units)
        
        columns = []
        
            
        for x in numbers:
            columns.append(map(str, x))
        
        rows = []
        for i in range(len(columns[0])):
            row = [x[i] for x in columns]
            rows.append(row)
            
        lines = map(lambda  x : self.column_separator.join(x), rows)
        
        for x in lines:
            self.stream.write(x)
            self.stream.write('\n')
            
        
    def write_footer(self):
        for x in self.footer_lines():
            self.stream.write(self.footer_prefix_string)
            self.stream.write(x)
            self.stream.write('\n')
        
    def header_lines(self):
        result = []
        result.append(' '.join(self.attribute_names))
        return result
        
    def footer_lines(self):
        return []
        
    def convert_string_to_number(self, string):
        return float(string)
        
    def new_set(self, number_of_items):
        return core.Particles(number_of_items)
        

        
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
        
            
        
    
