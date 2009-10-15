import re

from amuse.support.core import late, print_out

def strip_indent(string_with_indents):
    return re.sub('^ *\n', '', string_with_indents.rstrip())

class CreateDescriptionOfALegacyFunctionDefinition(object):

    @late
    def out(self):
        return print_out()
    
    def start(self):
        self.output_function_description()
        self.out.lf()
        self.output_cfunction_definition()
        self.output_fortran_function_definition()
        self.out.lf()
        self.output_parameter_descriptions()
        self.output_parameter_returntype()
        self.out.lf()

    
    def output_function_description(self):
        self.output_multiline_string(self.specification.description)
        self.out.lf()

    
    def output_multiline_string(self, string):
        lines = string.split('\n')
        first = True
        for line in lines:
            if first:
                first = False
            else:
                self.out.lf()
            self.out + line
            
    def output_cfunction_definition(self):
        self.out + '.. code-block:: c'
        self.out.indent()
        self.out.lf().lf()
        
        if self.specification.result_type is None:
            self.out + 'void '
        else:
            self.out + self.specification.result_type
            self.out + ' '
        self.out + self.specification.name
        self.out + '('
        self.out.indent()
        first = True
        for parameter in self.specification.parameters:
            if first:
                first = False
            else:
                self.out + ', '
            
            length_of_the_argument_statement = len(parameter.datatype) + len(parameter.name) + 3
            new_length_of_the_line = self.out.number_of_characters_on_current_line + length_of_the_argument_statement
            if new_length_of_the_line > 74:
                self.out.lf()
            self.out + parameter.datatype
            if parameter.is_output():
                self.out + ' *'
            self.out + ' ' + parameter.name
        self.out + ');'
        self.out.dedent()
        
        self.out.dedent()
        self.out.lf().lf()
        
    def output_fortran_function_definition(self):
        self.out + '.. code-block:: fortran'
        self.out.indent()
        self.out.lf().lf()
        
        if self.specification.result_type is None:
            type = 'SUBROUTINE'
        else:
            type = 'FUNCTION'

        self.out + type + ' '
        self.out + self.specification.name
        self.out + '('
        self.out.indent()
        self.out.indent()
        first = True
        for parameter in self.specification.parameters:
            if first:
                first = False
            else:
                self.out + ', '
            
            length_of_the_argument_statement = len(parameter.name)
            new_length_of_the_line = self.out.number_of_characters_on_current_line + length_of_the_argument_statement
            if new_length_of_the_line > 74:
                selt.out + ' &'
                self.out.lf()
            self.out + parameter.name
        
        self.out + ')'
        self.out.dedent()
        
        dtype_to_parameters = {}
        for parameter in self.specification.parameters: 
            parameters = dtype_to_parameters.get(parameter.datatype,[])
            parameters.append(parameter)
            dtype_to_parameters[parameter.datatype] = parameters
            
        dtype_to_fortan = {'int32':'INTEGER', 'float64':'DOUBLE PRECISION', 'float32':'REAL'}
        for dtype,parameters in dtype_to_parameters.iteritems():
            typestring = dtype_to_fortan[dtype]
            first = True
            
            self.out.lf()
            self.out + typestring + ' :: '
            
            for parameter in parameters:
                if first:
                    first = False
                else:
                    self.out + ', '
              
                length_of_the_argument_statement = len(parameter.name)
                new_length_of_the_line = self.out.number_of_characters_on_current_line + length_of_the_argument_statement
                if new_length_of_the_line > 74:
                    self.out.lf()
                    self.out + typestring + ' :: '
                
                self.out + parameter.name    
            
        if not self.specification.result_type is None:
            typestring = dtype_to_fortan[self.specification.result_type]
            self.out.lf()
            self.out + typestring + ' :: ' + self.specification.name

        self.out.dedent()
        self.out.lf()
        self.out + 'END ' + type
        
        self.out.dedent()
        self.out.lf().lf()
    
    def output_parameter_descriptions(self):
        for parameter in self.specification.parameters:
            self.out.lf()
            self.out + ':param ' + parameter.name + ': '
            self.out.indent()
            self.output_multiline_string(strip_indent(parameter.description))
            self.out.dedent()
            self.out.lf()
            self.out + ':type ' + parameter.name + ': '
            self.out + parameter.datatype + ', '
            self.output_parameter_direction(parameter) 


    def output_parameter_direction(self, parameter):
        #self.out + '('
        if parameter.direction == self.specification.IN:
            self.out + 'IN'
        if parameter.direction == self.specification.INOUT:
            self.out + 'INOUT'
        if parameter.direction == self.specification.OUT:
            self.out + 'OUT'
        #self.out + ')'

    def output_parameter_returntype(self):
        if self.specification.result_type is None:
            return
        self.out.lf()
        self.out + ':returns: '
        self.out.indent()
        self.output_multiline_string(strip_indent(self.specification.result_doc))
        self.out.dedent()
        
class CreateInterfaceDefinitionDocument(object):
    @late
    def out(self):
        return print_out()
        
    def start(self):
        pass
