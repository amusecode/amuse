from amuse.support.core import late
from amuse.support import exceptions
from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import dtypes
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.tools import create_definition
from amuse.rfi.core import LegacyFunctionSpecification

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('ints_in', 'ints_out',
                    'HEADER_INTEGER_COUNT', 'int', 'MPI_INT'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'HEADER_LONG_COUNT', 'long long int', 'MPI_LONG_LONG_INT'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'HEADER_FLOAT_COUNT', 'float', 'MPI_FLOAT'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'HEADER_DOUBLE_COUNT', 'double', 'MPI_DOUBLE'),
    'bool' : DTypeSpec('booleans_in', 'booleans_out',
                    'HEADER_BOOLEAN_COUNT', 'int', 'MPI_INTEGER'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'HEADER_STRING_COUNT', 'int', 'MPI_INTEGER'),
})

from amuse.support.core import late
from amuse.support import exceptions
from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import dtypes
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.tools.create_python_worker import CreateAPythonWorker
from amuse.rfi.tools import create_definition
from amuse.rfi.core import LegacyFunctionSpecification

import sys
import os
import inspect
class MakeCythonCodeString(GenerateASourcecodeString):
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
       
         


class GenerateACythonStringOfAFunctionSpecification(MakeCythonCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
   
        
    def output_function_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            #if not parameter.is_input():
            #    continue
                
            if first:
                first = False
            else:
                self.out + ', '
            name = 'in_' if parameter.name == 'in' else parameter.name                                                                                  
            if self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.IN:
                
                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                else:     
                    self.out + 'numpy.ndarray[' + spec.type + ', ndim=1, mode="c"] ' + name
            elif self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out + 'int ' + name
            else:
                self.out + name
                    

    def start(self):
        
        self.specification.prepare_output_parameters()
        
        self.out + 'def' + ' ' + self.specification.name + '('
        self.output_function_parameters()
        self.out + '):'

        
                    
        self.out.indent().lf()
        self.output_output_variables()
        #self.output_parameters()
        
        if self.specification.must_handle_array:
            first_input_parameter = self.first_input_parameter()
            if 0:
                length_parameter = self.length_parameter()
                self.out.lf() + 'cdef int LENGTH_' + length_parameter.name + ' = '                             
                if first_input_parameter.direction == LegacyFunctionSpecification.IN:
                    self.out + 'len('+first_input_parameter.name + ')'
                else:
                    self.out + 'len('+first_input_parameter.name + '.value)'
                self.out.lf()                                         
            
        self.output_function_start()
        
        self.output_function_call_parameters()
        self.output_function_end()
        self.output_filloutput_variables()
        #self.output_parameters()
        
        if 0:
            if self.specification.must_handle_array:
                if not self.specification.result_type is None:
                    spec = self.dtype_to_spec[self.specification.result_type]
                    self.out.lf() + 'for (int i = 1 ; i < call_count; i++){'
                    self.out.indent()
                    self.out.lf() + spec.output_var_name + '[i]' + ' = ' + spec.output_var_name + '[0]' + ';'
                    self.out.dedent()
                    self.out.lf() + '}'
            elif self.specification.can_handle_array:
                self.out.dedent()
                self.out.lf() + '}'
                
        if not self.specification.result_type is None:
            self.out.lf() + 'return __result__'
            
            
        self.out.dedent().lf()
        
        self._result = self.out.string
    


    def index_string(self, index, must_copy_in_to_out=False):
        if self.specification.must_handle_array and not must_copy_in_to_out:
            if index == 0:
                return '0'
            else:
                return '( %d * call_count)' % index
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                return '( %d * call_count) + i' % index
        else:
            return index
    
    
    def input_var(self, name, index):
        if self.specification.must_handle_array:
            self.output_var(name, index)
        else:
            self.out.n() + name
            self.out + '[' + self.index_string(index) + ']'
        
    def output_var(self, name, index):
        self.out.n() + '&' + name
        self.out + '[' + self.index_string(index) + ']'
    
    def output_output_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            
            if self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.INOUT:
                
                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                else:     
                    self.out.lf() + 'cdef numpy.ndarray[' + spec.type + ', ndim=1, mode="c"]  inout_' + parameter.name + ' = ' +  parameter.name + '.value'
            elif self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.OUT:
                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                
                spec = self.dtype_to_spec[parameter.datatype]
            
                
                self.out.lf() + 'cdef numpy.ndarray[' + spec.type + ', ndim=1, mode="c"]  output_' + parameter.name + ' = ' +  'numpy.zeros('+self.length_parameter().name +', dtype = '+self.numpy_dtype(spec)+')'
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                spec = self.dtype_to_spec[parameter.datatype]
            
                self.out.lf() + 'cdef '
                
                if parameter.datatype == 'string':
                    self.out + 'char *'
                else:
                    self.out + spec.type                                                                                         
                self.out + ' ' + 'output_' +  parameter.name
            elif parameter.direction == LegacyFunctionSpecification.INOUT:
                spec = self.dtype_to_spec[parameter.datatype]
            
                if parameter.datatype == 'string':
                    self.out.lf() + 'py_' + '_byte_' +  parameter.name  + ' = ' + parameter.name + '.value' + ".encode('UTF-8')"
                self.out.lf() + 'cdef '
                
                if parameter.datatype == 'string':
                    self.out + 'char *'
                else:
                    self.out + spec.type                                     
                self.out + ' ' + 'inout_' +  parameter.name + ' = '                                                 
                if parameter.datatype == 'string':
                    self.out + 'py_' + '_byte_' +  parameter.name
                else:
                    self.out + parameter.name + '.value'
                
        


    def output_copy_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'for (int i = 0 ; i < call_count; i++){'
                    self.out.indent()

                self.out.n() + spec.output_var_name
                self.out + '[' + self.index_string(parameter.output_index, must_copy_in_to_out=True) + ']'
                self.out + ' = '
                self.out + spec.input_var_name + '[' + self.index_string(parameter.input_index, must_copy_in_to_out=True) + ']' + ';'
            
                if self.specification.must_handle_array:
                    self.out.dedent()
                    self.out.lf() + '}'
                                
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for parameter in self.specification.output_parameters:
            count = dtype_to_count.get(parameter.datatype, 0)
            dtype_to_count[parameter.datatype] = count + 1
                
        if not self.specification.result_type is None:
            count = dtype_to_count.get(self.specification.result_type, 0)
            dtype_to_count[self.specification.result_type] = count + 1
            
        for dtype in dtype_to_count:       
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_count[dtype]
            self.out.n() 
            self.out + 'header_out[' + spec.counter_name 
            self.out + '] = ' + count + ' * call_count;'
            pass
            
    def output_function_end(self):
        #if len(self.specification.parameters) > 0:
        #    self.out.n()
            
        self.out + ')' + ';'
        

    def output_function_start(self):
        self.out.n()                                     
        
        

        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
             
            self.out + 'cdef '
            if self.specification.result_type == 'string':
                self.out + 'char *'
            else:
                self.out + spec.type                                     
            self.out + ' __result__'
            self.out + ' = '
        self.out + 'c_' + self.specification.name + '('
        


    def output_casestmt_start(self):
        self.out + 'case ' + self.specification.id + ':'
        
    def output_casestmt_end(self):
        self.out.n() + 'break;'
        
        

    def output_function_call_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ', '
            name = 'in_' if parameter.name == 'in' else parameter.name                                                        
            if self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.IN:

                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                else:     
                    self.out + '&' + name + '[0]'
            elif self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.OUT:

                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                else:     
                    self.out + '&output_' + name + '[0]'
            elif self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.INOUT:

                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                else:     
                    self.out + '&inout_' + name + '[0]'
            else:
                #if parameter.direction == LegacyFunctionSpecification.LENGTH:
                #     name = 'LENGTH_' + name
                if parameter.direction == LegacyFunctionSpecification.OUT:
                    self.out + '&output_'
                elif parameter.direction == LegacyFunctionSpecification.INOUT:
                    self.out + '&inout_'
                self.out + name
                if parameter.datatype == 'string' and parameter.direction == LegacyFunctionSpecification.IN:
                    self.out + ".encode('UTF-8')"
                    


    def output_parameters(self):
        for parameter in self.specification.parameters:
            name = 'in_' if parameter.name == 'in' else parameter.name
            self.out.lf() + 'print ' + name
        self.out.lf()
                
        

    def output_filloutput_variables(self):
        for parameter in self.specification.parameters:
            if parameter.direction == LegacyFunctionSpecification.OUT:
                self.out.lf() +  parameter.name + '.value = ' +  'output_' + parameter.name
                if parameter.datatype == 'string':
                    self.out + ".decode('UTF-8')"
            elif parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() +  parameter.name + '.value = ' +  'inout_' + parameter.name

                if parameter.datatype == 'string':
                    self.out + ".decode('UTF-8')"
                
                
        


    def first_input_parameter(self):
        
        for parameter in self.specification.parameters:
            if parameter.is_input():
                return parameter

    def length_parameter(self):
        
        for parameter in self.specification.parameters:
            if parameter.direction == LegacyFunctionSpecification.LENGTH:
                return parameter

    def numpy_dtype(self, spec):
        
        ctype = spec.type
        if ctype == 'int':
            return 'numpy.int32'
        elif ctype == 'long long int':
            return 'numpy.int64'
        else:
            return 'numpy.'+ctype

class GenerateACythonDefinitionStringFromAFunctionSpecification(MakeCythonCodeString):
   
        
    def start(self):
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self._result = self.out.string
            
    def output_function_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ', '
                
            if parameter.datatype == 'string':
                self.out + 'char'
            else:
                self.out + spec.type
            if parameter.is_output() or (parameter.is_input() and self.specification.must_handle_array):
                self.out + ' ' + '*'
                if parameter.datatype == 'string':
                    self.out  + '*'
            else:
                if parameter.datatype == 'string':
                    self.out  + ' ' + '*'
                
            
    def output_function_end(self):
        self.out + ')' + ';'
        
    def output_function_start(self):
        self.out.n()
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.type
            self.out + ' '
        else:
            self.out + 'void' + ' '
        self.out + 'c_'+self.specification.name +' "'+self.specification.name +'" '+ '('
        
class GenerateACythonSourcecodeStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def specification_class(self):
        raise exceptions.AmuseException("No specification_class set, please set the specification_class first")
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def output_sourcecode_for_function(self):
        return GenerateACythonStringOfAFunctionSpecification()
    
    def start(self):
        self.out + 'import numpy'
        self.out.lf() + 'cimport numpy'

        self.output_local_includes()

        self.output_mpi_defs()

        self.out.lf() + 'cdef extern from "worker_code.h":'
        
        self.out.indent().lf()
        
        self.output_definitions_for_functions()
        
        self.out.lf().dedent()
        
        self.output_sourcecode_for_functions()
        
        #self.out + self.new_executable_script_string()
        
        self._result = self.out.string
        


    def output_definitions_for_functions(self):
        for x in self.interface_functions:
            if x.specification.id == 0:
                continue
            if not self.must_include_interface_function_in_output(x):
                continue
                
            self.out.lf()
            uc = GenerateACythonDefinitionStringFromAFunctionSpecification()
            uc.specification = x.specification
            uc.out = self.out
            uc.start()
            self.out.lf()
    
    
    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        #for cls in self.ignore_functions_from_specification_classes:
        #    if hasattr(cls, x.specification.name):
        #        return False
        
        return True
    @late
    def template_string(self):
        path = self.template_dir
        path = os.path.join(path, 'cython_code_script.template')
            
        with open(path, "r") as f:
            template_string = f.read()
        
        return template_string

    @late
    def template_dir(self):
        return os.path.dirname(__file__)

    def new_executable_script_string(self):
        return self.template_string.format(
            syspath = ','.join(map(repr, sys.path)),
            worker_module = 'interface',
            interface_module = inspect.getmodule(self.specification_class).__name__,
            interface = self.specification_class.__name__,
        )



    def output_local_includes(self):
        if hasattr(self.specification_class, 'include_headers'):
            for x in self.specification_class.include_headers:
                self.out.n() + 'cdef extern from "' + x + '":'
                self.out.n() + '    '+'pass'
        self.out.lf()

    def output_mpi_defs(self):
        self.out.lf() + "cimport mpi4py.MPI"
        self.out.lf() + 'cdef extern from "mpi.h":\n    pass'
        self.out.lf() + 'cdef extern from "amuse_mpi.h":'
        self.out.lf() + "    " + 'int c_set_comm_world "set_comm_world" (mpi4py.MPI.MPI_Comm world)'
        self.out.lf().lf() + "def set_comm_world(mpi4py.MPI.Comm comm not None):"
        self.out.lf() + "    " + 'return c_set_comm_world(comm.ob_mpi)'
        self.out.lf()


class GenerateACHeaderStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def ignore_functions_from_specification_classes(self):
        return []
        
    @late
    def underscore_functions_from_specification_classes(self):
        return []
        
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
        
    @late
    def make_extern_c(self):
        return True
    
    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        for cls in self.ignore_functions_from_specification_classes:
            if hasattr(cls, x.specification.name):
                return False
        
        return True
        
    def output_sourcecode_for_function(self):
        return GenerateACHeaderDefinitionStringFromAFunctionSpecification()
        
    def start(self):  
        if self.make_extern_c:
            self.out + 'extern "C" {'
            self.out.indent().lf()
            
        self.output_sourcecode_for_functions()
        
        if self.make_extern_c:
            self.out.dedent().lf() + '}'
        
        self.out.lf()
        
        
        self._result = self.out.string
        

class GenerateACStubStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec
        
    @late
    def make_extern_c(self):
        return False
    
    def output_sourcecode_for_function(self):
        return create_definition.CreateCStub()

    def must_include_interface_function_in_output(self, x):
        return not x.specification.name.startswith("internal__")
     
    def start(self):  
    
        self.output_local_includes()
        
        self.out.lf()
        
        if self.make_extern_c:
            self.out + 'extern "C" {'
            self.out.indent().lf()
            
        self.output_sourcecode_for_functions()
        
        if self.make_extern_c:
            self.out.dedent().lf() + '}'
        
        self.out.lf()
        
        self._result = self.out.string
        
    
    def output_local_includes(self):
        self.out.n()
        if hasattr(self.specification_class, 'include_headers'):
            for x in self.specification_class.include_headers:
                self.out.n() + '#include "' + x + '"'
    
        
        
        
class GenerateACythonStartScriptStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def specification_class(self):
        raise exceptions.AmuseException("No specification_class set, please set the specification_class first")
    
   
    def start(self):
        self.out + self.new_executable_script_string()
        
        self._result = self.out.string
        
    
    @late
    def template_string(self):
        path = self.template_dir
        path = os.path.join(path, 'cython_code_script.template')
            
        with open(path, "r") as f:
            template_string = f.read()
        
        return template_string

    @late
    def template_dir(self):
        return os.path.dirname(__file__)

    def new_executable_script_string(self):
        return self.template_string.format(
            executable = sys.executable,
            syspath = ','.join(map(repr, sys.path)),
            worker_module = self.worker_module,
            interface_module = inspect.getmodule(self.specification_class).__name__,
            interface = self.specification_class.__name__,
        )





    @late
    def worker_module(self):
        return self.name_of_outputfile



    @late
    def name_of_outputfile(self):
        return 'interface'



