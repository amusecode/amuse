"""
"""
from amuse.support.core import late
from amuse.support.core import print_out

from zlib import crc32

class LegacyCall(object):
    """A legacy_call implements the runtime call to the remote process.
    """
    def __init__(self, interface, owner, specification):
        self.interface = interface
        self.owner = owner
        self.specification = specification
    
    def __call__(self, *arguments_list, **keyword_arguments):
        dtype_to_values_and_keyword = {
            'd' : ([],'doubles_in'),
            'i' : ([],'ints_in')
        }
        names_in_argument_list = set([])
        for index, argument in argument_list:
            name, dtype, direction = self.specification.parameters[index]
            names_in_argument_list.add(name)
            values = dtype_to_values_and_keyword[dtype]
            values.append(argument)
            
            
        call_keyword_arguments = {}
        for values, keyword in dtype_to_values_and_keyword.values():
            call_keyword_arguments[keyword] = values
            
        (ints_out, doubles_out) = self.do_call(self.specification.id, **call_keyword_arguments)
        return (ints_out, doubles_out)

class legacy_function(object):
    """The meta information for a function call to a code
    """
    def __init__(self, specification_function):
        self.specification_function = specification_function
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return self.new_legacy_call(instance, owner)
        
    def __set__(self, instance, value):
        return
        
    def legacy_call(self, instance, owner):
        return None
        
    def new_legacy_call(self, instance, owner):
        return LegacyCall(instance, owner, self.specification)
    
    def to_c_string(self):
        uc = MakeACStringOfALegacyFunctionSpecification()
        uc.specification = self.specification
        return uc.result 
        
    @late
    def specification(self):
        result = self.specification_function()
        if result.name is None:
            result.name = self.specification_function.__name__
        if result.id is None:
            result.id = crc32(result.name)
        return result
        
class RemoteFunction(object):
    IN = object()
    OUT = object()
    INOUT = object()
    
    def __init__(self):
        self.parameters = []
        self.name = None
        self.id = None
        self.result_type = None
        
    def addParameter(self, name, dtype = 'i', direction = IN):
        self.parameters.append((name, dtype, direction))
    
    
        
                
        
        
class MakeACStringOfALegacyFunctionSpecification(object):
    def __init__(self):
        pass
    
    @late  
    def result(self):
        self.start()
        return self._result
    
    @late
    def dtype_to_spec(self):
        class DTypeSpec(object):
            def __init__(self, c_input_var_name, c_output_var_name, c_output_counter_name):
                self.number_of_inputs = 0
                self.number_of_outputs = 0
                self.c_input_var_name = c_input_var_name
                self.c_output_var_name = c_output_var_name
                self.c_output_counter_name = c_output_counter_name
             
        return {
            'i' : DTypeSpec('ints_in','ints_out', 'number_of_ints'),
            'd' : DTypeSpec('doubles_in', 'doubles_out','number_of_doubles')}
        
    def start(self):
        self.out = print_out()
        self.output_casestmt_start()
        self.out.indent()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        self.output_lines_with_number_of_outputs()
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
        
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for name, dtype, direction in self.specification.parameters:
            spec = self.dtype_to_spec[dtype]
            
            if first:
                first = False
            else:
                self.out + ' ,'
                
            if direction == RemoteFunction.IN:
                self.out.n() + spec.c_input_var_name + '[' + spec.number_of_inputs + ']'
                spec.number_of_inputs += 1
            if direction == RemoteFunction.INOUT:
                self.out.n() + '&' + spec.c_input_var_name + '[' + spec.number_of_inputs + ']'
                spec.number_of_inputs += 1
            elif direction == RemoteFunction.OUT:
                self.out.n() + '&' + spec.c_output_var_name + '[' + spec.number_of_outputs + ']'
                spec.number_of_outputs += 1
    
        self.out.dedent()
    def output_lines_with_inout_variables(self):
        dtype_to_incount = {}
        
        for name, dtype, direction in self.specification.parameters:
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_incount.get(dtype, 0)
        
            if direction == RemoteFunction.IN:
                dtype_to_incount[dtype] = count + 1
            if direction == RemoteFunction.INOUT:
                self.out.n() + spec.c_output_var_name + '[' + spec.number_of_outputs + ']' + ' = ' + spec.c_input_var_name + '[' + count + ']'+';'
                spec.number_of_outputs += 1
                dtype_to_incount[dtype] = count + 1
    
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT or direction == RemoteFunction.INOUT:
                count = dtype_to_count.get(dtype, 0)
                dtype_to_count[dtype] = count + 1
                
        if not self.specification.result_type is None:
            count = dtype_to_count.get(dtype, 0)
            dtype_to_count[dtype] = count + 1
            
        for dtype in dtype_to_count:       
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_count[dtype]
            self.out.n() 
            self.out + 'reply.' + spec.c_output_counter_name + ' = ' + count + ';'
            pass
    def output_function_end(self):
        self.out.n() + ')' + ';'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.c_output_var_name+ '[' + spec.number_of_outputs + ']' + ' = '
            spec.number_of_outputs += 1
        self.out + self.specification.name + '('
    def output_casestmt_start(self):
        self.out + 'case ' + self.specification.id + ':'
    def output_casestmt_end(self):
        self.out.n() + 'break;'
        

        
            