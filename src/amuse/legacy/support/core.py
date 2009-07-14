"""
"""
from amuse.support.core import late
from amuse.support.core import print_out

class legacy_call(object):
    """A legacy_call implements the runtime call to the remote process.
    """
    def __init__(self, interface):
        self.interface = interface
    
    def __call__(self, *arguments_list, **keyword_arguments):
        pass

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
        return legacy_call(instance)
    
    def to_c_string(self):
        uc = MakeACStringOfALegacyFunctionSpecification()
        uc.specification = self.specification
        return uc.result 
        
    @late
    def specification(self):
        return self.specification_function()
        
class RemoteFunction(object):
    IN = object()
    OUT = object()
    INOUT = object()
    
    def __init__(self):
        self.parameters = []
        self.name = "<noname>"
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
            elif direction == RemoteFunction.OUT:
                self.out.n() + '&' + spec.c_output_var_name + '[' + spec.number_of_outputs + ']'
                spec.number_of_outputs += 1
    
        self.out.dedent()
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT:
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
        

        
            