from amuse.support.core import late
from amuse.support import exceptions
import numpy

from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import dtypes
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.tools import create_definition
from amuse.rfi.core import LegacyFunctionSpecification
dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('integers_in','integers_out', 
                    'number_of_integers', 'integer', 'MPI_INTEGER'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'number_of_doubles', 'real*8', 'MPI_REAL8'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'number_of_floats', 'real*4', 'MPI_REAL4'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'number_of_strings', 'integer', 'MPI_INTEGER'),
    'bool' : DTypeSpec('bools_in', 'bools_out',
                    'number_of_bools', 'logical', 'MPI_LOGICAL'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'number_of_longs', 'integer*8', 'MPI_LONG_LONG_INT'),
})
        
redirect_outputs_function_template = """
function internal__redirect_outputs(stdoutfile, stderrfile)
    implicit none
    
    character(LEN=*) , INTENT(IN) :: stdoutfile, stderrfile
    character(1024) :: fullname
    integer :: mpi_rank1, mpi_err1, internal__redirect_outputs
    
    CLOSE(UNIT=5) ! always close stdin
    
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank1, mpi_err1)
    
    if (stdoutfile .NE. 'none' ) then
        
        close(UNIT=6)
        if (stdoutfile .NE. '/dev/null') then
            write (fullname, '(A,".",I3.3)')  stdoutfile, mpi_rank1
            open(unit=6, file=trim(fullname), access="append")
        end if
    end if
    
    if (stderrfile .NE. 'none') then
        close(UNIT=0)
        
        if (stderrfile .NE. '/dev/null') then
            write( fullname, '(A,".",I3.3)' )  stderrfile, mpi_rank1
            open(unit=0, file=trim(fullname), access="APPEND")
        end if
        
    end if
    
    internal__redirect_outputs = 0
end function
"""

class GenerateAFortranStringOfAFunctionSpecification(GenerateASourcecodeString):
    MAX_STRING_LEN = 256
    
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
    
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
        
    def index_string(self, index, must_copy_in_to_out = False):
        if self.specification.must_handle_array and not must_copy_in_to_out:
            if index == 0:
                return '1'
            else:
                return '( %d * len_in) + 1' % (index )
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                if index == -1:
                    return "i - 1"
                else:
                    return '( %d * len_in) + i' % index
        else:
            return index + 1
            
    def start(self):        
        self.specification.prepare_output_parameters()
         
        self.output_casestmt_start()
        self.out.indent()
        
        self.output_lines_before_with_clear_out_variables()
        self.output_lines_before_with_clear_input_variables()
        
        if self.specification.must_handle_array:
            pass
        elif self.specification.can_handle_array:
            self.out.lf() + 'do i = 1, len_in, 1'
            self.out.indent()
        
        self.output_lines_before_with_inout_variables()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'DO i = 2, len_in'
                self.out.indent()
                self.out.lf() + spec.output_var_name + '(i)' + ' = ' + spec.output_var_name + '(1)'
                self.out.dedent()
                self.out.lf() + 'END DO'
        elif self.specification.can_handle_array:
            self.out.dedent()
            self.out.lf() + 'end do'
            
        self.output_lines_with_number_of_outputs()
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
        
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
                self.out + ' &'
            else:
                self.out + ' ,&'
                
            if parameter.direction == LegacyFunctionSpecification.IN:
                if parameter.datatype == 'string':
                    self.out.n() + 'input_characters('
                    self.out  + '( (' + self.index_string(parameter.input_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.input_index) + ')* ' + self.MAX_STRING_LEN + ') +'
                    self.out  +  '(' + spec.input_var_name + '(' + self.index_string(parameter.input_index) + ')' + '-' 
                    self.out  + 'get_offset(' + self.index_string(parameter.input_index) + ' - 1 , '+spec.input_var_name +') ))'
                    self.out  + ')'
                else:
                    self.out.n() + spec.input_var_name 
                    self.out + '(' + self.index_string(parameter.input_index) + ')'
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string':
                    self.out.n() + 'output_characters('
                    self.out  + '((' + self.index_string(parameter.output_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.output_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
                    self.out  + ')'
                else:
                    self.out.n() + spec.input_var_name 
                    self.out + '(' + self.index_string(parameter.input_index) + ')'
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                if parameter.datatype == 'string':
                    self.out.n() + 'output_characters('
                    self.out  + '((' + self.index_string(parameter.output_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.output_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
                    self.out  + ')'
                else:
                    self.out.n() + spec.output_var_name
                    self.out + '(' + self.index_string(parameter.output_index) + ')'
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out.n() + 'len_in'
                
        self.out.dedent()
        
    def output_lines_with_inout_variables(self):
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'DO i = 1, len_in'
                    self.out.indent() 
                    
                self.out.n() + spec.output_var_name 
                self.out + '(' + self.index_string(parameter.output_index, must_copy_in_to_out = True)  + ')' 
                self.out + ' = ' 
                self.out + spec.input_var_name + '(' + self.index_string(parameter.input_index, must_copy_in_to_out = True) + ')'
        
                if self.specification.must_handle_array:
                    self.out.dedent() 
                    self.out.lf() + 'END DO'
    
    def output_lines_before_with_clear_out_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.is_output():
                if parameter.datatype == 'string': 
                    self.out.lf() + 'output_characters = "x"'  
                    return
     
    def output_lines_before_with_clear_input_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.is_input():
                if parameter.datatype == 'string': 
                    self.out.lf() + 'input_characters = "x"'  
                    return
     
                
                    
    def output_lines_before_with_inout_variables(self):
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            
            if parameter.direction == LegacyFunctionSpecification.IN:
                if parameter.datatype == 'string':
                    self.out.n() + 'input_characters('
                    self.out  + '( (' + self.index_string(parameter.input_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.input_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
                    self.out  + ') = &'
                    self.out.lf()
                    self.out + 'characters('
                    self.out + 'get_offset(' + self.index_string(parameter.input_index) + ' - 1 , '+spec.input_var_name +')'
                    self.out  + ':' + spec.input_var_name + '(' + self.index_string(parameter.input_index) + ')'
                    self.out  + ')' 
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string':
                    self.out.n() + 'output_characters('
                    self.out  + '( (' + self.index_string(parameter.output_index) + ')* ' + self.MAX_STRING_LEN + ')'
                    self.out  + ':' + '(((' + self.index_string(parameter.output_index) + ')+1) * ' + self.MAX_STRING_LEN + ' - 1)'
                    self.out  + ') = &'
                    self.out.lf()
                    self.out + 'characters('
                    self.out + 'get_offset(' + self.index_string(parameter.input_index) + ' - 1 , '+spec.input_var_name +')'
                    self.out  + ':' + spec.input_var_name + '(' + self.index_string(parameter.input_index) + ')'
                    self.out  + ')' 
                    
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
            self.out + spec.counter_name + '_out'
            self.out + ' = ' + count 
            pass
            
    def output_function_end(self):
        self.out + ' &'
        self.out.n() + ')'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            if self.specification.result_type == 'string':
                self.out + 'output_characters('
                self.out  + '( (' + self.index_string(0) + ')* ' + self.MAX_STRING_LEN + ')'
                self.out  + ':' + '(((' + self.index_string(0) + ')+1)*' + self.MAX_STRING_LEN + '-1)'
                self.out  + ') = &'
                self.out.lf()
            else:
                self.out + spec.output_var_name
                self.out + '(' + self.index_string(0) + ')' + ' = '
        else:    
            self.out + 'CALL ' 
        self.out +  self.specification.name + '('
        
    def output_casestmt_start(self):
        self.out + 'CASE(' + self.specification.id + ')'
        
    def output_casestmt_end(self):
        self.out.n() 
        
        
class MakeAFortranStringOfALegacyGlobalSpecification(GenerateASourcecodeString):
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
            
    def start(self):
        self.output_casestmt_start()
        self.out.indent()
        
        spec = self.dtype_to_spec[self.legacy_global.datatype]
        self.out.n() + 'if (' +spec.counter_name +'_in'
        self.out + ' == 1 ) then'
        self.out.indent()
        self.out.n() + self.legacy_global.name + ' = ' 
        self.out + spec.input_var_name  + '[1]'
        self.out.dedent()
        self.out.n() + 'else'
        self.out.indent()
        self.out.n() + spec.counter_name + '_out'
        self.out + ' = ' + 1 
        self.out.n() + spec.output_var_name + '[1]' 
        self.out + ' = ' + self.legacy_global.name
        self.out.dedent()
        self.out.n()
        
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
        
        
    def output_casestmt_start(self):
        self.out + 'CASE(' + self.legacy_global.id  + ')'
        
    def output_casestmt_end(self):
        self.out.n() 


class GenerateAFortranSourcecodeStringFromASpecificationClass(GenerateASourcecodeStringFromASpecificationClass):
    MAX_STRING_LEN = 256

    @late
    def dtype_to_spec(self):
        return dtype_to_spec 
   
    @late
    def number_of_types(self):
        return len(self.dtype_to_spec)
        
    @late
    def length_of_the_header(self):
        return 2 + self.number_of_types
        
        
    def output_sourcecode_for_function(self):
        return GenerateAFortranStringOfAFunctionSpecification()
   
    def start(self):
        self.output_character_index_function()
        self.output_runloop_function_def_start()
        self.output_switch_start()
        self.output_sourcecode_for_functions()
        self.output_switch_end() 
        
        
        self.output_runloop_function_def_end()
        self.output_main()
        self._result = self.out.string

    def output_mpi_include(self):
        self.out.n() + "INCLUDE 'mpif.h'"
        
  
            
    def output_modules(self):
        self.out.n()
        if hasattr(self.specification_class, 'use_modules'):
            for x in self.specification_class.use_modules:
                self.out.n() + 'use ' + x 
                
    def must_include_declaration_of_function(self, x):
        if x.specification.name.startswith("internal__"):
            return False
        
        return True
        
        
    def output_declarations_for_the_functions(self):
        if not hasattr(self.specification_class, 'use_modules'):
            for x in self.interface_functions:
                if not self.must_include_declaration_of_function(x):
                    continue
                    
                specification = x.specification
                if specification.id == 0:
                    continue
                if specification.result_type is None:
                    continue
                if specification.result_type == 'string':
                    type = 'CHARACTER(len=255)'
                else:
                    spec = self.dtype_to_spec[specification.result_type]
                    type = spec.type
                self.out.lf() +  type + ' :: ' + specification.name
        
    def output_allocate_arrays(self):
        maximum_number_of_inputvariables_of_a_type = 255
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'ALLOCATE('
                self.out + dtype_spec.input_var_name 
                self.out + '( maxlen *' + max + ')'
                self.out + ')'
            
            max =self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'ALLOCATE(' 
                self.out + dtype_spec.output_var_name 
                self.out + '( maxlen * ' + max + ')'
                self.out + ')'

        

    def output_character_index_function(self):
        self.out.lf().lf() + 'FUNCTION get_offset(index, offsets)'
        self.out.indent()
        self.out.lf() + 'integer :: index, get_offset'
        self.out.lf() + 'integer,dimension(:) :: offsets'
        self.out.lf() + 'IF (index .lt. 1) THEN'
        self.out.indent().lf() + 'get_offset = 1'
        self.out.dedent().lf() + 'ELSE'
        self.out.indent().lf() + 'get_offset = offsets(index) + 2'
        self.out.dedent().lf() + 'END IF'
        
        self.out.dedent()
        self.out.lf() + 'END FUNCTION'
        
    def output_character_index_inderface(self):
        self.out.lf() + 'INTERFACE'
        self.out.indent()
        self.out.lf().lf() + 'FUNCTION get_offset(index, offsets)'
        self.out.indent()
        self.out.lf() + 'integer :: index, get_offset'
        self.out.lf() + 'integer,dimension(:) :: offsets'
        self.out.dedent().lf() + 'END FUNCTION'
        self.out.dedent().lf() + 'end interface'
        
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'SUBROUTINE run_loop'
        self.out.indent()
        self.output_modules()
        self.out.lf().lf() + 'IMPLICIT NONE'
        self.out.lf()
        self.output_mpi_include()
        self.out.n() + 'integer :: rank, parent, ioerror, maxlen = 255'
        self.out.n() + 'integer :: must_run_loop'
        self.out.n() + 'integer i, str_len, offset'
        
        self.out.n() + 'character (len=100000) :: characters'
        self.out.n() + 'character (len=100000) :: output_characters'
        self.out.n() + 'character (len=100000) :: input_characters'
        
        #self.out.n() + 'integer mpiStatus(MPI_STATUS_SIZE,4)'
        self.out.lf().lf() + 'integer header('
        self.out + self.length_of_the_header + ')'
        self.out.lf().lf() + 'integer :: tag_in, tag_out'
        self.out.lf().lf() + 'integer :: len_in, len_out'
        self.output_character_index_inderface()
        self.out.lf()
        self.output_declarations_for_the_functions()
        
        for dtype_spec in self.dtype_to_spec.values():
            self.out.lf() + dtype_spec.type
            self.out + ', DIMENSION(:), ALLOCATABLE :: '
            self.out + dtype_spec.input_var_name
            
            self.out.lf() + dtype_spec.type
            self.out + ', DIMENSION(:), ALLOCATABLE :: '
            self.out + dtype_spec.output_var_name
            
            self.out.lf() + 'integer ::' + ' '
            self.out + dtype_spec.counter_name + '_out'
            self.out + ', ' + dtype_spec.counter_name + '_in'
            
        
        self.out.lf()
        self.output_allocate_arrays()
        self.out.lf()
        self.out.lf().lf() + 'call MPI_COMM_GET_PARENT(parent, ioerror)'
        self.out.lf()      + 'call MPI_COMM_RANK(parent, rank, ioerror)'
        self.out.lf().lf() + 'must_run_loop = 1'
        self.out.lf().lf() + 'do while (must_run_loop .eq. 1)'
        self.out.indent()
       
       
        self.out.lf() + 'call MPI_BCast(header, '
        self.out + self.length_of_the_header
        self.out + ', MPI_INTEGER, 0,'
        self.out + ' parent,&'
        self.out.indent().lf() + 'ioerror)'
        self.out.dedent()
        self.out.lf().lf() + 'tag_in = header(1)'
        self.out.lf().lf() + 'len_in = header(2)'
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + spec.counter_name + '_in' + ' =  '
            self.out + 'header(' + (i+3) + ')'
        
        self.out.lf().lf() + 'tag_out = tag_in'
        self.out.lf() + 'len_out = len_in'
        
          
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + spec.counter_name + '_out' + ' =  0'
        self.out.lf()
        
        
        self.out.lf() + 'IF (len_in .gt. maxlen) THEN'
        self.out.indent()
        self.out.lf() + 'maxlen = len_in + 255;'
        self.output_deallocate_statements()
        self.output_allocate_arrays()
        self.out.dedent()
        self.out.lf() + 'END IF'
    
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'if (' + spec.counter_name + '_in'
            self.out + ' .gt. 0) then'
    
            self.out.indent().lf() + 'call MPI_BCast('
            self.out + spec.input_var_name + ', '
            self.out + spec.counter_name + '_in'
            self.out + ' * ' + 'len_in'
            self.out + ', &'
    
            self.out.indent().n() + spec.mpi_type
            self.out + ', 0, parent,&'
            self.out.n() + 'ioError);'
            if dtype == 'string':
            
                self.out.dedent()
                
                self.out.lf() + 'call MPI_BCast('
                self.out + 'characters' + ', '
                self.out + spec.input_var_name + '('+spec.counter_name + '_in' + '* len_in' +') + 1'
                self.out + ', &'
    
                self.out.indent().n() + 'MPI_CHARACTER'
                self.out + ', 0, parent,&'
                self.out.n() + 'ioError);'
                
            self.out.dedent().dedent().lf()
            self.out + 'end if'
         
    def output_switch_start(self):
        self.out.lf().lf() + 'SELECT CASE (tag_in)'
        self.out.indent()
        self.out.lf() + 'CASE(0)'
        self.out.indent().lf()+'must_run_loop = 0'
        self.out.dedent()
        
    def output_switch_end(self):
        self.out.lf() + 'CASE DEFAULT'
        self.out.indent().lf() + 'tag_out = -1'
        self.out.dedent()
        self.out.dedent().lf() + 'END SELECT'
        
    def output_runloop_function_def_end(self):
        
        self.out.lf().lf() + 'if (rank .eq. 0 ) then'
        self.out.indent().lf()
        
        self.out.lf().lf() + 'header(1) = tag_out'
        self.out.lf() + 'header(2) = len_out'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf()  + 'header(' + (i+3) + ')'
            self.out + ' = ' + spec.counter_name + '_out'
        
        self.out.lf().lf() + 'call MPI_SEND(header, '
        self.out + self.length_of_the_header
        self.out + ', MPI_INTEGER,'
        self.out + ' 0, 999, &'
        self.out.indent().lf() + 'parent, ioerror);'
        self.out.dedent().lf()
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            
            self.out.lf() + 'if (' + spec.counter_name + '_out'
            self.out + ' .gt. 0) then'
            self.out.indent().lf()
            if dtype == 'string':
                self.out.lf() + 'offset = 1'
                self.out.lf() + 'DO i = 1, '+spec.counter_name + '_out * len_out' + ',1'
                self.out.indent().lf()
                self.out.lf() + 'str_len = LEN_TRIM(output_characters(i*' + self.MAX_STRING_LEN + ':((i+1)*' + self.MAX_STRING_LEN + ')-1))'
                self.out.lf() + 'output_characters(offset:offset+255) = output_characters(i*' + self.MAX_STRING_LEN + ':((i+1)*' + self.MAX_STRING_LEN + ')-1)'
                self.out.lf() + 'offset = offset + str_len - 1'
                self.out.lf() + spec.output_var_name + '(i) = offset'
                self.out.lf() + 'offset = offset + 2'
                self.out.dedent().lf() + 'END DO'
                
                self.out.lf() + 'call MPI_SEND('
                self.out + spec.output_var_name + ', ' +  spec.counter_name + '_out'
                self.out + ' * len_out'
                self.out + ', &'
                self.out.indent().lf() + spec.mpi_type + ', 0, 999, &'
                self.out.lf() + 'parent, ioerror)'
                self.out.dedent()
                
                self.out.lf() + 'call MPI_SEND('
                self.out + 'output_characters' + ', ' + 'offset -1'
                self.out + ', &'
                self.out.indent().lf() + 'MPI_CHARACTER' + ', 0, 999, &'
                self.out.lf() + 'parent, ioerror)'
                
                
            else:
                self.out.lf()
                self.out + 'call MPI_SEND('
                self.out + spec.output_var_name + ', ' +  spec.counter_name + '_out'
                self.out + ' * len_out'
                self.out + ', &'
                self.out.indent().lf() + spec.mpi_type + ', 0, 999, &'
                self.out.lf() + 'parent, ioerror)'
                
            self.out.dedent().dedent().lf() +'end if'
        
        
        
        self.out.dedent().lf()+ 'end if'
        self.out.indent().lf().lf()
            
        self.out.dedent()
        self.out.lf() + 'end do'
        
        
        self.out.lf()
        self.output_deallocate_statements()
            
        self.out.lf() + 'call MPI_COMM_DISCONNECT(parent, ioerror)'
        self.out.lf() + 'return'
        
        
        self.out.lf().lf() + 'CONTAINS'
        self.out.indent()
        self.output_redirect_output()
        self.out.dedent()
        
        self.out.dedent().lf() + 'end subroutine'
    
    def output_deallocate_statements(self):
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'DEALLOCATE(' + dtype_spec.input_var_name  + ')'
            
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'DEALLOCATE(' + dtype_spec.output_var_name  + ')'
    
    def output_redirect_output(self):
        self.out.lf().lf() + redirect_outputs_function_template
        
    def output_main(self):
        self.out.lf().lf() + 'program muse_worker'
        self.out.indent()
        self.output_mpi_include()
        self.out.lf() + 'integer :: provided,ioerror'
        self.out.lf() + 'call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ioerror)'
        self.out.lf().lf() + 'call run_loop()'
        self.out.lf().lf() + 'call MPI_FINALIZE(ioerror)'
        self.out.dedent().lf()+'end program muse_worker'
        
        
        

class GenerateAFortranStubStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec
  
    @late
    def ignore_functions_from_specification_class(self):
        return []
        
    def output_sourcecode_for_function(self):
        result = create_definition.CreateFortranStub()
        result.output_definition_only = False
        return result
        
    def start(self):  
    
        self.output_modules()
        
        self.out.lf()
        
        self.output_sourcecode_for_functions()
        
        self.out.lf()
        
        self._result = self.out.string
        
    
    def must_include_interface_function_in_output(self, x):
        if x.specification.name.startswith("internal__"):
            return False
            
        for cls in self.ignore_functions_from_specification_class:
            if hasattr(cls, x.specification.name):
                return False
        
        return True
        
    def output_modules(self):
        self.out.n()
        if hasattr(self.specification_class, 'use_modules'):
            for x in self.specification_class.use_modules:
                self.out.n() + 'use ' + x 
        
    
        
        

        
       
    
        
        
        
