import re

from amuse.support.core import late, print_out

def strip_indent(string_with_indents):
    return re.sub('^ *\n', '', string_with_indents.rstrip())

class CreateDescriptionOfALegacyFunctionDefinition(object):

    @late
    def out(self):
        return print_out()
    
    def start(self):
        self.output_cfunction_definition()
        self.out.indent()
        self.out.lf().lf()
        self.out.lf().lf()
        self.output_function_description()
        self.output_parameter_descriptions()
        self.output_parameter_returntype()
        self.out.dedent()
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
        self.out + '.. function:: '
        if self.specification.result_type is None:
            self.out + 'void '
        else:
            self.out + self.specification.result_type
            self.out + ' '
        self.out + self.specification.name
        self.out + '('
        first = True
        for parameter in self.specification.parameters:
            if first:
                first = False
            else:
                self.out + ', '
            self.out + parameter.datatype
            if parameter.is_output():
                self.out + ' *'
            self.out + ' ' + parameter.name
        self.out + ')'
    
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
