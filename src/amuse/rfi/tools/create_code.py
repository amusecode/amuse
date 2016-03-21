import numpy
import sys
import os

from amuse.support.core import late, print_out
from amuse.rfi.core import legacy_function

class DTypeSpec(object):
    def __init__(self, input_var_name, output_var_name, counter_name, 
        type, mpi_type = 'UNKNOWN'):
        self.input_var_name = input_var_name
        self.output_var_name = output_var_name
        self.counter_name = counter_name
        self.type = type
        self.mpi_type = mpi_type
        
    


dtypes = ['int32', 'int64', 'float32', 'float64', 'bool', 'string']

class GenerateASourcecodeString(object):
    _result = None

    def __init__(self):
        pass
    
    @late  
    def result(self):
        if self._result is None:
            self.start()
        return self._result
        
    @late
    def out(self):
        return print_out()
        
    @late
    def must_generate_mpi(self):
        if 'CFLAGS' in os.environ:
            return not (os.environ['CFLAGS'].find('-DNOMPI') >= 0)
        else:
            return True


class GenerateASourcecodeStringFromASpecificationClass(GenerateASourcecodeString):
    
    @late
    def interface_functions(self):
        attribute_names = dir(self.specification_class)
        interface_functions = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(self.specification_class, x)
            if isinstance(value, legacy_function):
                interface_functions.append(value)
        interface_functions.sort(key= lambda x: x.specification.id)
        return interface_functions
        
    @late
    def mapping_from_dtype_to_maximum_number_of_inputvariables(self):
        result = None
        for x in self.interface_functions:
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
        for x in self.interface_functions:
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
                
    def must_include_interface_function_in_output(self, x):
        return True
        
    def output_sourcecode_for_functions(self):
        for x in self.interface_functions:
            if x.specification.id == 0:
                continue
            if not self.must_include_interface_function_in_output(x):
                continue
                
            self.out.lf()
            uc = self.output_sourcecode_for_function()
            uc.specification = x.specification
            uc.out = self.out
            uc.start()
            self.out.lf()
    
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
        
        
