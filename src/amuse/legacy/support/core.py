"""
"""
from amuse.support.core import late
from amuse.support.core import print_out

class legacy_call(object):
    def __init__(self, interface):
        self.interface = interface
    
    def __call__(self, *arguments_list, **keyword_arguments):
        pass

class legacy_function(object):
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
class MakeACStringOfALegacyFunctionSpecification(object):
    def __init__(self):
        pass
    
    @late  
    def result(self):
        self.start()
        return self._result
    
    @late
    def dtype_to_spec(self):
        return {
            'i' : [0,'ints_in', 0 , 'ints_out'],
            'd' : [0,'doubles_in', 0 , 'doubles_out'],
        }
    def start(self):
        self.out = print_out()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
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
                self.out.n() + spec[1] + '[' + spec[0] + ']'
                spec[0] += 1
            elif direction == RemoteFunction.OUT:
                self.out.n() + '&' + spec[3] + '[' + spec[2] + ']'
                spec[2] += 1
    
        self.out.dedent()
    def output_function_end(self):
        self.out.n() + ')' + ';'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec[3] + '[' + spec[2] + ']' + ' = '
            spec[2] += 1
        self.out + self.specification.name + '('
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
    
    
        
                
        
        