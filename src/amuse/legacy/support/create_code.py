from amuse.support.core import late, print_out
from amuse.legacy.support.core import legacy_function, legacy_global

class DTypeSpec(object):
    def __init__(self, input_var_name, output_var_name, output_counter_name, type):
        self.number_of_inputs = 0
        self.number_of_outputs = 0
        self.input_var_name = input_var_name
        self.output_var_name = output_var_name
        self.output_counter_name = output_counter_name
        self.type = type

class MakeCodeString(object):
    def __init__(self):
        pass
    
    @late  
    def result(self):
        self.start()
        return self._result
        
    @late
    def out(self):
        return print_out()

class MakeCodeStringOfAClassWithLegacyFunctions(MakeCodeString):
    
    @late
    def legacy_functions(self):
        attribute_names = dir(self.class_with_legacy_functions)
        legacy_functions = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(self.class_with_legacy_functions, x)
            if isinstance(value, legacy_function):
                legacy_functions.append(value)
        
        legacy_functions.sort(key= lambda x: x.specification.id)
        return legacy_functions
        
    @late
    def legacy_globals(self):
        attribute_names = dir(self.class_with_legacy_functions)
        result = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(self.class_with_legacy_functions, x)
            if isinstance(value, legacy_global):
                result.append(value)
        
        result.sort(key= lambda x: x.id)
        return result
           
    def output_legacy_functions(self):
        for x in self.legacy_functions:
            if x.specification.id == 0:
                continue
            self.out.lf()
            uc = self.make_legacy_function()
            uc.specification = x.specification
            uc.out = self.out
            uc.start()
    
    def output_legacy_globals(self):
        for x in self.legacy_globals:
            self.out.lf()
            uc = self.make_legacy_global()
            uc.legacy_global = x
            uc.out = self.out
            uc.start()
            
    def output_extra_content(self):
        self.out.lf()
        if hasattr(self.class_with_legacy_functions, 'extra_content'):
            self.out.n() + self.class_with_legacy_functions.extra_content
    
            
    
    
        
        
        
