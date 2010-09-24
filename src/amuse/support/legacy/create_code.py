from amuse.support.core import late, print_out
from amuse.support.legacy.core import legacy_function, legacy_global

import numpy
import sys
class DTypeSpec(object):
    def __init__(self, input_var_name, output_var_name, counter_name, 
        type, mpi_type = 'UNKNOWN'):
        self.input_var_name = input_var_name
        self.output_var_name = output_var_name
        self.counter_name = counter_name
        self.type = type
        self.mpi_type = mpi_type

dtypes = ['float64', 'int32', 'float32', 'string', 'bool']

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
            if isinstance(value, legacy_function) and self.handle_legacy_function(x, value):
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
        
    
    @late
    def mapping_from_dtype_to_maximum_number_of_inputvariables(self):
        result = None
        for x in self.legacy_functions:
            local = {}
            for parameter in x.specification.input_parameters:
                count = local.get(parameter.datatype, 0)
                local[parameter.datatype] = count + 1
            
            
            if result is None:
                result = local
            else:
                for key, count in local.iteritems():
                    previous_count = result.get(key, 0)
                    result[key] = max(count, previous_count)
                    
        return result
                
    @late
    def mapping_from_dtype_to_maximum_number_of_outputvariables(self):
        result = None
        for x in self.legacy_functions:
            local = {}
            for parameter in x.specification.output_parameters:
                count = local.get(parameter.datatype, 0)
                local[parameter.datatype] = count + 1
                
            if not x.specification.result_type is None:
                count = local.get(x.specification.result_type, 0)
                local[x.specification.result_type] = count + 1
            
            if result is None:
                result = local
            else:
                for key, count in local.iteritems():
                    previous_count = result.get(key, 0)
                    result[key] = max(count, previous_count)
                    
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
            self.out.lf()
    
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
    
    def handle_legacy_function(self, name_of_the_function, function):
        return True
    
    
class DTypeToSpecDictionary(object):
    
    def __init__(self, dict):
        self.mapping = {}
        for datatype, value in dict.iteritems():
            self.mapping[datatype] = value
        
    def __getitem__(self, datatype):
        return self.mapping[datatype]
    
    def __len__(self):
        return len(self.mapping)
        
    def values(self):
        return self.mapping.values()
        
    def keys(self):
        return self.mapping.keys()
        
        
