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
        x = CreateFortranStub()
        x.out = self.out
        x.specification = self.specification
        x.start()
        
    
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
        
class CreateFortranStub(object):
    @late
    def out(self):
        return print_out()
        
    def start(self):
        self.out + '.. code-block:: fortran'
        self.out.indent()
        self.out.lf().lf()
        self.output_subprogram_start()
        self.output_parameter_type_definiton_lines()
        self.output_subprogram_end()
                           
        self.out.dedent()
        self.out.lf().lf()

    
    @late 
    def specification_is_for_function(self):
        return not self.specification.result_type is None
    
    @late 
    def subprogram_string(self):
        if self.specification_is_for_function:
            return 'FUNCTION'
        else:
            return 'SUBROUTINE'
            
    @late
    def dtype_to_parameters(self):
        result = {}
        for parameter in self.specification.parameters: 
            parameters = result.get(parameter.datatype,[])
            parameters.append(parameter)
            result[parameter.datatype] = parameters
        return result
        
    @late
    def dtype_to_fortantype(self):
        return {
            'int32':'INTEGER' , 
            'float64':'DOUBLE PRECISION' , 
            'float32':'REAL' ,
            'string':'CHARACTER(LEN=*)',
        }
        
    def output_subprogram_start(self):
        self.out + self.subprogram_string + ' '
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
                self.out + ' &'
                self.out.lf()
            self.out + parameter.name
        
        self.out + ')'
        self.out.dedent()

    def output_parameter_type_definiton_lines(self):
        for dtype,parameters in self.dtype_to_parameters.iteritems():
            typestring = self.dtype_to_fortantype[dtype]
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
    
    def output_function_type(self):
        if self.specification_is_for_function:
            typestring = self.dtype_to_fortantype[self.specification.result_type]
            self.out.lf()
            self.out + typestring + ' :: ' + self.specification.name

    def output_subprogram_end(self):
        self.out.dedent()
        self.out.lf()
        self.out + 'END ' + self.subprogram_string

        