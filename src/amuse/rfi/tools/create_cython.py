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
                    'HEADER_BOOLEAN_COUNT', 'bool', 'MPI_C_BOOL'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'HEADER_STRING_COUNT', 'char *', 'MPI_INTEGER'),
})

INTERFACE_DEFINITION = """
    interface
        type(C_ptr) &
        function C_malloc(size) bind(C,name="malloc")
          import c_size_t, c_ptr
          integer(C_size_t), value, intent(in) :: size
        end function C_malloc
    end interface
"""

STRING_UTILITY_FUNCTIONS = """
    
    subroutine C_F_string_ptr(C_string, F_string)
        type(C_ptr), intent(in) :: C_string
        character(len=*), intent(out) :: F_string
        character(len=1,kind=C_char), dimension(:), pointer :: p_chars
        integer :: i
        F_string = ' '
        if (.not. C_associated(C_string)) then
          F_string = ' '
        else
          call C_F_pointer(C_string,p_chars,[huge(0)])
          i=1
          do while(p_chars(i)/=NUL .and. i<=len(F_string))
            F_string(i:i) = p_chars(i)
            i=i+1
          end do
          if (i<len(F_string)) F_string(i:) = ' '
        end if
    end subroutine C_F_string_ptr

    subroutine F_C_string_ptr(F_string, C_string, C_string_len)
        character(len=*), intent(in) :: F_string
        type(C_ptr), intent(in) :: C_string ! target = intent(out)
        integer, intent(in), optional :: C_string_len  ! Max string length,
                                                       ! INCLUDING THE TERMINAL NUL
        character(len=1,kind=C_char), dimension(:), pointer :: p_chars
        integer :: strlen
        strlen = len(F_string)
        if (present(C_string_len)) then
          if (C_string_len <= 0) return
          strlen = min(strlen,C_string_len-1)
        end if
        if (.not. C_associated(C_string)) then
          return
        end if
        call C_F_pointer(C_string,p_chars,[strlen+1])
    !   do i=1,strlen
          p_chars(1)(1:strlen) = F_string(1:strlen)
    !   enddo
        p_chars(strlen+1) = NUL
    end subroutine F_C_string_ptr
"""

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


    @late
    def dtype_to_fortran_type(self):
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
                    self.out + 'numpy.ndarray[' + spec.type + ', ndim=1, mode="c", cast=True] ' + name
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
                    self.out.lf() + 'cdef numpy.ndarray[' + spec.type + ', ndim=1, mode="c", cast=True]  inout_' + parameter.name + ' = ' +  parameter.name + '.value'
            elif self.specification.must_handle_array and parameter.direction == LegacyFunctionSpecification.OUT:
                if parameter.datatype == 'string':
                    raise Exception("unknown...")
                
                spec = self.dtype_to_spec[parameter.datatype]
            
                
                self.out.lf() + 'cdef numpy.ndarray[' + spec.type + ', ndim=1, mode="c", cast=True]  output_' + parameter.name + ' = ' +  'numpy.zeros('+self.length_parameter().name +', dtype = '+self.numpy_dtype(spec)+')'
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
                    self.out.lf() + 'py_' + '_byte_' +  parameter.name  + ' = ' + parameter.name + '.value' #+ ".encode('UTF-8')"
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
                #if parameter.datatype == 'string' and parameter.direction == LegacyFunctionSpecification.IN:
                #    self.out + ".encode('UTF-8')"
                    


    def output_parameters(self):
        for parameter in self.specification.parameters:
            name = 'in_' if parameter.name == 'in' else parameter.name
            self.out.lf() + 'print ' + name
        self.out.lf()
                
        

    def output_filloutput_variables(self):
        for parameter in self.specification.parameters:
            if parameter.direction == LegacyFunctionSpecification.OUT:
                self.out.lf() +  parameter.name + '.value = ' +  'output_' + parameter.name
                #if parameter.datatype == 'string':
                #    self.out + ".decode('UTF-8')"
            elif parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() +  parameter.name + '.value = ' +  'inout_' + parameter.name

                #if parameter.datatype == 'string':
                #    self.out + ".decode('UTF-8')"
                
                
        


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

        self.out + 'c_'+self.specification.name +' "' + self.function_name_prefix + self.specification.name +'" '+ '('
        



    @late
    def function_name_prefix(self):
        return ''



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
        self.out.lf() + 'cdef extern from "stdbool.h":'
        self.out.lf() + '  ctypedef bint bool' 
        self.output_local_includes()

        self.output_mpi_defs()
        
        if self.include_filename_for_functions is None:
            self.out.lf() + 'cdef extern:'
        else:
            self.out.lf() + 'cdef extern from "{0}":'.format(self.include_filename_for_functions)
        
        self.out.indent().lf()
        
        self.output_definitions_for_functions()
        
        self.out.lf().dedent()
        
        self.output_sourcecode_for_functions()
        
        #self.out + self.new_executable_script_string()
        
        self._result = self.out.string
        




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


    @late
    def include_filename_for_functions(self):
        if hasattr(self.specification_class, 'include_headers') and len(self.specification_class.include_headers) >= 1:
            return self.specification_class.include_headers[0]
        return None



    def output_definitions_for_functions(self):
        for x in self.interface_functions:
            if x.specification.id == 0:
                continue
            if not self.must_include_interface_function_in_output(x):
                continue
                
            self.out.lf()
            uc = GenerateACythonDefinitionStringFromAFunctionSpecification()
            uc.specification = x.specification
            uc.function_name_prefix = self.function_name_prefix
            uc.out = self.out
            uc.start()
            self.out.lf()

    @late
    def function_name_prefix(self):
        return ''


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
        return self.cython_import if self.cython_import else self.name_of_outputfile




    @late
    def name_of_outputfile(self):
        return 'interface'



    @late
    def cython_import(self):
        return ""


class GenerateAFortranInterfaceStringOfAFunctionSpecification(MakeCythonCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
   
        
    def output_function_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ', &\n &'
            name = parameter.name                                                                                      
            self.out + name
                    


    def start(self):
        
        self.specification.prepare_output_parameters()
        
        self.output_function_start()
        self.output_function_parameters()
        self.out + ') &\n &'
        
        if not self.specification.result_type is None: 
            self.out + ' result(rrreeesss) &\n &'
            
        self.out + ' ' + 'bind(c, name = "' + self.function_name_prefix + self.specification.name + '")'
        self.out.n()
        self.out.indent().lf()
        self.out.lf() + 'implicit none'
        self.output_function_parameter_definitions()
        self.output_function_string_definitions()
        self.output_function_logical_definitions()
        if self.specification.result_type == 'string':
            self.out.lf() + 'character(len=4096) :: string_rrreeesss'
        if self.has_string_output_parameters:
            self.out.lf() + 'integer(8) :: sz'
            
        
        if not self.has_modules and not self.specification.result_type is None:
            fortran_type = self.dtype_to_fortran_simple_type[self.specification.result_type]
            self.out.lf() + fortran_type + ' :: '
            self.out + self.specification.name
        
        self.output_function_logical_input_copies()
        self.output_function_string_input_copies()
        self.output_function_string_output_clears()
        self.output_function_call_start()
        self.output_function_call_parameters()
        self.output_function_call_end()
        self.output_function_string_output_copies()
        self.output_function_logical_output_copies()
       
        self.out.dedent().lf()
        self.output_function_end()
        
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
    
    def output_function_end(self):
       
        if not self.specification.result_type is None:
            self.out.lf() + 'END FUNCTION ' + self.function_name
        else:
            self.out.lf() + 'END SUBROUTINE ' + self.function_name
            



    def output_function_start(self):
        self.out.n()                                       

        if not self.specification.result_type is None:
            fortran_type = self.dtype_to_fortran_type[self.specification.result_type]      
            self.out + fortran_type
                            
            self.out +  ' function '  + self.function_name + '('
        else:
            self.out +  ' subroutine '  + self.function_name + '('
        





    def output_function_call_parameters(self):        
        first = True
        
        for parameter in self.specification.parameters:
            if first:
                first = False
            else:
                self.out + ', &\n &'
            name = parameter.name                          
            if parameter.datatype == 'string':
                self.out + 'string_' + name              
            elif parameter.datatype == 'bool':
                self.out + 'logical_' + name
            else:     
                self.out + name
                    





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
            
            




    @late
    def dtype_to_fortran_type(self):
        return  {
            'int32' : 'INTEGER(kind = c_int)',
            'int64' : 'INTEGER(kind = c_long)',
            'float32' : 'REAL(kind=c_float)',
            'float64' : 'REAL(kind=c_double)',
            'bool' : 'LOGICAL(kind = c_bool)',
            'string' : 'type(C_ptr)'
        }



    def output_function_parameter_definitions(self):        
        for parameter in self.specification.parameters:
            fortran_type = self.dtype_to_fortran_type[parameter.datatype]
            self.out.lf() + fortran_type                                                                     
            if  parameter.direction == LegacyFunctionSpecification.IN:
                self.out + ', intent(in)'
                if not self.specification.must_handle_array:
                    self.out + ', value'
            if  parameter.direction == LegacyFunctionSpecification.OUT:
                self.out + ', intent(out)'
            if  parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out + ', intent(inout)'
            if  parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out + ', intent(in), value'
            
            if self.specification.must_handle_array and not parameter.direction == LegacyFunctionSpecification.LENGTH:
                if parameter.datatype == 'string':
                    raise Exception("must_handle_array and string arguments not supported")
                self.out + ', dimension('+self.length_parameter().name+')'
            self.out + ' :: ' + parameter.name
                    





    @late
    def has_modules(self):
        return False
    def output_function_call_start(self):
        self.out.n()                                           

        if self.specification.result_type == 'string':
            self.out + 'string_rrreeesss = '
        elif not self.specification.result_type is None:
            self.out + 'rrreeesss = '
        else:
            self.out + 'call '
                            
        self.out + self.specification.name + '('
        





    def output_function_call_end(self):
        self.out + ')'

    @late
    def dtype_to_fortran_simple_type(self):
        return  {
            'int32' : 'INTEGER',
            'int64' : 'LONG',
            'float32' : 'REAL*4',
            'float64' : 'REAL*8',
            'bool' : 'LOGICAL',
            'string' : 'CHARACTER(len=4096)'
        }



    def output_function_string_definitions(self):        
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'string':
                continue
            if  parameter.direction == LegacyFunctionSpecification.IN:
                self.out.lf() + 'character(len=4096) :: string_'+parameter.name
            if  parameter.direction == LegacyFunctionSpecification.OUT:
                self.out.lf() + 'character(len=4096) :: string_'+parameter.name
            if  parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() + 'character(len=4096) :: string_'+parameter.name
            




    @late
    def has_string_parameters(self):
        for parameter in self.specification.parameters:
            if parameter.datatype == 'string':
                return True
        return False
    def output_function_string_input_copies(self):        
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'string':
                continue
            if  parameter.direction == LegacyFunctionSpecification.IN or parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() + 'call C_F_string_ptr({0}, string_{0})'.format(parameter.name)
                #self.out.lf() + 'print *, "::", trim(string_{0})'.format(parameter.name)
                pass
            





    def output_function_string_output_copies(self): 
        index = 1                                       
        if self.specification.result_type == 'string':             
            self.out.lf() + 'if (.NOT. C_associated(string_buffer{0})) then'.format(index)
            self.out.lf() + '   sz = 4097'
            self.out.lf() + '   string_buffer{0} = C_malloc(sz)'.format(index)
            self.out.lf() + 'end if'
            #self.out.lf() + 'print * , trim(string_rrreeesss)'
            self.out.lf() + 'call F_C_string_ptr(trim(string_{0}), string_buffer{1}, 4096)'.format('rrreeesss', index)
            self.out.lf() + 'rrreeesss = ' + 'string_buffer{0}'.format(index)
            index = index + 1                                                                 
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'string':
                continue
            if  parameter.direction == LegacyFunctionSpecification.OUT or parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() + 'if (.NOT. C_associated(string_buffer{0})) then'.format(index)
                self.out.lf() + '   sz = 4097'
                self.out.lf() + '   string_buffer{0} = C_malloc(sz)'.format(index)
                self.out.lf() + 'end if'
                #self.out.lf() + 'print * , trim(string_{0})'.format(parameter.name)
                self.out.lf() + 'call F_C_string_ptr(trim(string_{0}), string_buffer{1}, 4096)'.format(parameter.name, index)
                self.out.lf() + parameter.name + ' = ' + 'string_buffer{0}'.format(index)
                index = index + 1
            





    @late
    def crc32(self):
        try:
        
            from zlib import crc32
            try:
                if crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return crc32
            except Exception:
                #python 3, crc32 needs bytes...
                def python3_crc32(x):
                    x = crc32(bytes(x, 'ascii'))
                    return x - ((x & 0x80000000) <<1)
                if python3_crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return python3_crc32
        except Exception:
            pass
        try:
            from binascii import crc32
            try:
                if crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return crc32
            except Exception:
                #python 3, crc32 needs bytes...
                def python3_crc32(x):
                    x = crc32(bytes(x, 'ascii'))
                    return x - ((x & 0x80000000) <<1)
                if python3_crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return python3_crc32
        except Exception:
            pass
        
        raise Exception("No working crc32 implementation found!")


    @late
    def function_name(self):
        return "c_" + hex(abs(self.crc32(self.specification.name)))[2:]

    @late
    def has_string_output_parameters(self):
        if self.specification.result_type == 'string':
            return True                  
        for parameter in self.specification.parameters:
            if parameter.datatype == 'string' and (parameter.direction == LegacyFunctionSpecification.OUT or parameter.direction == LegacyFunctionSpecification.INOUT):
                return True
        return False


    def output_function_string_output_clears(self): 
        index = 1                                           
        if self.specification.result_type == 'string':             
            self.out.lf() + "string_{0} = ' '".format('rrreeesss')
                                                            
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'string':
                continue
            if  parameter.direction == LegacyFunctionSpecification.OUT:
                self.out.lf() + "string_{0} = ' '".format(parameter.name)
            




    def output_function_logical_definitions(self):        
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'bool':
                continue
            if self.specification.must_handle_array:
                self.out.lf() + 'LOGICAL, dimension('+self.length_parameter().name+') :: logical_'+parameter.name
            else:
                self.out.lf() + 'LOGICAL :: logical_'+parameter.name
            




    def output_function_logical_input_copies(self):        
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'bool':
                continue
            if  parameter.direction == LegacyFunctionSpecification.IN or parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() + 'logical_{0} = {0}'.format(parameter.name)
                




    def output_function_logical_output_copies(self):        
        for parameter in self.specification.parameters:
            if not parameter.datatype == 'bool':
                continue
            if  parameter.direction == LegacyFunctionSpecification.OUT or parameter.direction == LegacyFunctionSpecification.INOUT:
                self.out.lf() +' {0} = logical_{0}'.format(parameter.name)



    @late
    def function_name_prefix(self):
        return "ci_"



class GenerateAFortranInterfaceSourcecodeStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def specification_class(self):
        raise exceptions.AmuseException("No specification_class set, please set the specification_class first")
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def start(self):
        self.out + 'module CInterfaceModule'
        self.out.lf() + 'use :: iso_c_binding'
        for i in range(self.number_of_string_output_variables):
            self.out.lf() + 'type(C_ptr) string_buffer{0}'.format(i+1)
        if self.has_functions_with_strings:
            self.out.lf() + 'character(len=1,kind=C_char), parameter :: NUL = C_NULL_char'

        self.output_modules()
        
        if self.has_functions_with_strings:
            self.output_interface_definition()

        self.out.lf() + 'contains'

        
        if self.has_functions_with_strings:
            self.output_string_utility_functions()

        self.output_definitions_for_functions()
        
        self.out + 'end module CInterfaceModule'

        self._result = self.out.string
        







    def output_definitions_for_functions(self):
        for x in self.interface_functions:
            if x.specification.id == 0:
                continue
            if not self.must_include_interface_function_in_output(x):
                continue
                
            self.out.lf()
            uc = GenerateAFortranInterfaceStringOfAFunctionSpecification()
            uc.specification = x.specification
            uc.out = self.out
            uc.has_modules = self.has_modules
            uc.function_name_prefix = self.function_name_prefix
            uc.start()
            self.out.lf()
    
    


    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        #for cls in self.ignore_functions_from_specification_classes:
        #    if hasattr(cls, x.specification.name):
        #        return False
        
        return True
    def output_modules(self):
        self.out.n()
        if hasattr(self.specification_class, 'use_modules'):
            for x in self.specification_class.use_modules:
                self.out.n() + 'use ' + x     
    @late
    def has_modules(self):
        return hasattr(self.specification_class, 'use_modules') and len(self.specification_class.use_modules) > 0

    def output_string_utility_functions(self):
        self.out.n()
        self.out + STRING_UTILITY_FUNCTIONS
        self.out.n()

    @late
    def has_functions_with_strings(self):
        return self.mapping_from_dtype_to_maximum_number_of_outputvariables.get('string', 0) > 0 or self.mapping_from_dtype_to_maximum_number_of_inputvariables.get('string', 0) > 0         



    @late
    def number_of_string_output_variables(self):
        return self.mapping_from_dtype_to_maximum_number_of_outputvariables.get('string', 0) + 1


    def output_interface_definition(self):
        self.out.n()
        self.out + INTERFACE_DEFINITION
        self.out.n()
    @late
    def function_name_prefix(self):
        return "ci_"


class GenerateACFFIStringOfAFunctionSpecification(MakeCythonCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")


    @late
    def dtype_to_fortran_type(self):
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
                self.out + 'output_' + parameter.name + ' = '
                if parameter.datatype == 'string':
                    todo     #self.out + 'char *'
                else:
                    self.out +  'ffi.new("'+spec.type+' *")'
            elif parameter.direction == LegacyFunctionSpecification.INOUT:
                spec = self.dtype_to_spec[parameter.datatype]
            
                if parameter.datatype == 'string':
                    self.out.lf() + 'py_' + '_byte_' +  parameter.name  + ' = ' + parameter.name + '.value' #+ ".encode('UTF-8')"
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
                #if parameter.datatype == 'string' and parameter.direction == LegacyFunctionSpecification.IN:
                #    self.out + ".encode('UTF-8')"
                    


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
            
            





class GenerateAPythonStubStringOfAFunctionSpecification(MakeCythonCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")


    @late
    def dtype_to_fortran_type(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
   
        
    def output_function_parameters(self):        
        self.out + 'self'
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            #if not parameter.is_input():
            #    continue
                
            self.out + ', '
            name = parameter.name                                                                                          
            self.out + name
                    


    def start(self):
        
        self.specification.prepare_output_parameters()
        
        self.out + 'def' + ' ' + self.specification.name + '('
        self.output_function_parameters()
        self.out + '):'

        self.out.indent().lf()

        for parameter in self.specification.parameters:
            if parameter.is_output():
                self.out.lf() + parameter.name + '.value = ' + 0
        
                
        if not self.specification.result_type is None:
            self.out.lf() + 'return 0'
            
            
        self.out.dedent().lf()
        
        self._result = self.out.string
    



