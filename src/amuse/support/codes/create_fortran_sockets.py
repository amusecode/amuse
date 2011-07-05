from amuse.support.core import late
from amuse.support import exceptions
from amuse.support.codes.core import LegacyFunctionSpecification
from amuse.support.codes.create_code import GenerateASourcecodeString
from amuse.support.codes.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.support.codes.create_code import DTypeSpec, DTypeToSpecDictionary
from amuse.support.codes import create_definition

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('integers_in','integers_out', 
                    'HEADER_INTEGER_COUNT', 'integer (c_int32_t)', 'integer'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'HEADER_LONG_COUNT', 'integer (c_int64_t)', 'long'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'HEADER_FLOAT_COUNT', 'real (c_float)', 'float'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'HEADER_DOUBLE_COUNT', 'real (c_double)', 'double'),
    'bool' : DTypeSpec('booleans_in', 'booleans_out',
                    'HEADER_BOOLEAN_COUNT', 'logical (c_bool)', 'boolean'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'HEADER_STRING_COUNT', 'integer (c_int32_t)', 'integer'),
})

dtypes = ['int32', 'int64', 'float32', 'float64', 'bool', 'string']


forsockets_module_code = """module forsockets
    integer HEADER_FLAGS, HEADER_CALL_ID, HEADER_FUNCTION_ID, HEADER_CALL_COUNT, & 
        HEADER_INTEGER_COUNT, HEADER_LONG_COUNT, HEADER_FLOAT_COUNT, & 
        HEADER_DOUBLE_COUNT, HEADER_BOOLEAN_COUNT, HEADER_STRING_COUNT, & 
        HEADER_SIZE

    parameter (HEADER_FLAGS=1, HEADER_CALL_ID=2, HEADER_FUNCTION_ID=3, & 
        HEADER_CALL_COUNT=4, HEADER_INTEGER_COUNT=5, HEADER_LONG_COUNT=6, & 
        HEADER_FLOAT_COUNT=7, HEADER_DOUBLE_COUNT=8, & 
        HEADER_BOOLEAN_COUNT=9, HEADER_STRING_COUNT=10, & 
        HEADER_SIZE=10)

    interface
        subroutine receive_integers & 
            (ints, length) & 
            bind(c, name='forsockets_receive_integers')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: ints
            integer (c_int32_t), value :: length
        end subroutine receive_integers

        subroutine receive_longs & 
            (longs, length) & 
            bind(c, name='forsockets_receive_longs')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: longs
            integer (c_int32_t), value :: length
        end subroutine receive_longs

        subroutine receive_floats & 
            (floats, length) & 
            bind(c, name='forsockets_receive_floats')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: floats
            integer (c_int32_t), value :: length
        end subroutine receive_floats

        subroutine receive_doubles & 
            (doubles, length) & 
            bind(c, name='forsockets_receive_doubles')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: doubles
            integer (c_int32_t), value :: length
        end subroutine receive_doubles

        subroutine receive_booleans & 
            (booleans, length) & 
            bind(c, name='forsockets_receive_booleans')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: booleans
            integer (c_int32_t), value :: length
        end subroutine receive_booleans

        subroutine receive_string & 
            (string, length) & 
            bind(c, name='forsockets_receive_string')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: string
            integer (c_int32_t), value :: length
        end subroutine receive_string

        subroutine send_integers & 
            (ints, length) & 
            bind(c, name='forsockets_send_integers')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: ints
            integer (c_int32_t), value :: length
        end subroutine send_integers

        subroutine send_longs & 
            (longs, length) & 
            bind(c, name='forsockets_send_longs')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: longs
            integer (c_int32_t), value :: length
        end subroutine send_longs

        subroutine send_floats & 
            (floats, length) & 
            bind(c, name='forsockets_send_floats')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: floats
            integer (c_int32_t), value :: length
        end subroutine send_floats

        subroutine send_doubles & 
            (doubles, length) & 
            bind(c, name='forsockets_send_doubles')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: doubles
            integer (c_int32_t), value :: length
        end subroutine send_doubles

        subroutine send_booleans & 
            (booleans, length) & 
            bind(c, name='forsockets_send_booleans')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: booleans
            integer (c_int32_t), value :: length
        end subroutine send_booleans

        subroutine send_string & 
            (string, length) & 
            bind(c, name='forsockets_send_string')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: string
            integer (c_int32_t), value :: length
        end subroutine send_string

        subroutine forsockets_init & 
            (port) & 
            bind(c, name='forsockets_init')
            use iso_c_binding
            implicit none
            integer (c_int32_t), value :: port
        end subroutine forsockets_init

        subroutine forsockets_close & 
            () & 
            bind(c, name='forsockets_close')
            use iso_c_binding
            implicit none
        end subroutine forsockets_close

    end interface
end module forsockets
"""

redirect_outputs_function_code = """
function internal__redirect_outputs(stdoutfile, stderrfile)
    use iso_c_binding
    
    implicit none
    
    character(kind=c_char, len = *) , intent(in) :: stdoutfile, stderrfile
    integer(c_int32_t) :: internal__redirect_outputs
    
    print*, 'NOT redirecting output to', stdoutfile, ' and ', stderrfile
    call flush()
    
    internal__redirect_outputs = 0
end function
"""

redirect_outputs_interface_code = """
  interface
    integer (c_int32_t) function internal__redirect_outputs(stdoutfile, stderrfile)
        use iso_c_binding
        character(kind=c_char, len = *) , intent(in) :: stdoutfile, stderrfile
    end function
  end interface
"""

main_program_code = """
  program amuse_worker
    use iso_c_binding
    use forsockets
    
    implicit none
    
    include 'mpif.h'
    integer :: provided,ioerror, port
    character(len=32) :: port_string

    call mpi_init_thread(mpi_thread_multiple, provided, ioerror)

    call get_command_argument(1, port_string)

    read (port_string,*) port

    call forsockets_init(port)
    
    call run_loop()
    
    call mpi_finalize(ioerror)

    call forsockets_close()

  end program amuse_worker
"""

                
string_receive_code = """
      call receive_integers(c_loc(strings_in), header_in(HEADER_STRING_COUNT))

      !print*, 'received string header:', strings_in
      !call flush()


      do i = 1, header_in(HEADER_STRING_COUNT), 1
        if (strings_in(i) .gt. MAX_STRING_LENGTH) then
            print*, 'error! cannot receive strings exeeding length ', MAX_STRING_LENGTH
        end if
      end do

      !space for all strings in this one call
      allocate(characters_in(header_in(HEADER_STRING_COUNT)))

      do i = 1, header_in(HEADER_STRING_COUNT), 1
          characters_in(i) = ' '
          
          call receive_string(c_loc(characters_in(i)), strings_in(i))

          print*, 'received string: ', characters_in(i)

          call flush()

      end do
"""

string_send_code = """
      !figure out length of all strings
      do i = 1, header_out(HEADER_STRING_COUNT), 1
          strings_out(i) = len_trim(characters_out(i))
      end do

      !send string header
      call send_integers(c_loc(strings_out), header_out(HEADER_STRING_COUNT))

      do i = 1, header_out(HEADER_STRING_COUNT), 1
          print*, 'sending string ', characters_out(i)
          call flush()
          call send_string(c_loc(characters_out(i)), strings_out(i))
      end do     
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
                return '( %d * length) + 1' % (index )
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                if index == -1:
                    return "i - 1"
                else:
                    return '( %d * length) + i' % index
        else:
            return index + 1
            
    def start(self):        
        self.specification.prepare_output_parameters()
         
        self.output_casestmt_start()
        self.out.indent()

        self.output_lines_with_number_of_outputs()
        
        self.output_lines_before_with_clear_out_variables()
#        self.output_lines_before_with_clear_input_variables()

        if self.specification.must_handle_array:
            pass
        elif self.specification.can_handle_array:
            self.out.lf() + 'do i = 1, length, 1'
            self.out.indent()
            
        
#        self.output_lines_before_with_inout_variables()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'DO i = 2, length'
                self.out.indent()
                self.out.lf() + spec.output_var_name + '(i)' + ' = ' + spec.output_var_name + '(1)'
                self.out.dedent()
                self.out.lf() + 'END DO'
        elif self.specification.can_handle_array:
            self.out.dedent()
            self.out.lf() + 'end do'
            

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
                    self.out.n() + 'trim(characters_in'
                    self.out + '(' + self.index_string(parameter.input_index) + '))'      
                else:
                    self.out.n() + spec.input_var_name 
                    self.out + '(' + self.index_string(parameter.input_index) + ')'
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string':
                    self.out.n() + 'trim(characters_in'
                    self.out + '(' + self.index_string(parameter.input_index) + '))'      
                else:
                    self.out.n() + spec.input_var_name 
                    self.out + '(' + self.index_string(parameter.input_index) + ')'
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                if parameter.datatype == 'string':
                    self.out.n() + 'characters_out'
                    self.out + '(' + self.index_string(parameter.output_index) + ')'      
                else:
                    self.out.n() + spec.output_var_name
                    self.out + '(' + self.index_string(parameter.output_index) + ')'
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out.n() + 'length'
                
        self.out.dedent()
        
    def output_lines_with_inout_variables(self):
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'DO i = 1, length'
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
                    self.out.lf() + 'allocate(characters_out(header_out(HEADER_STRING_COUNT)))'
                    self.out.lf() + "characters_out = ' '"  
  
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
            self.out + 'header_out(' + spec.counter_name + ')'
            self.out + ' = ' + count + " * length" 
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
        self.output_forsockets_module()
        self.output_redirect_output()
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
                
    def output_declarations_for_the_functions(self):
        if not hasattr(self.specification_class, 'use_modules'):
            for x in self.interface_functions:
                specification = x.specification
                if specification.id == 0:
                    continue
                if specification.result_type is None:
                    continue
                if specification.name == 'internal__redirect_outputs':
                    continue
                if specification.result_type == 'string':
                    type = 'CHARACTER(len=255)'
                else:
                    spec = self.dtype_to_spec[specification.result_type]
                    type = spec.type
                self.out.lf() +  type + ' :: ' + specification.name
        
    def output_allocate_arrays(self):
        maximum_number_of_inputvariables_of_a_type = 255
        
        for i, dtype in enumerate(dtypes):
            dtype_spec = self.dtype_to_spec[dtype]
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'allocate('
                self.out + dtype_spec.input_var_name 
                self.out + '( max_length *' + max + ')'
                self.out + ')'
            
            max =self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'allocate(' 
                self.out + dtype_spec.output_var_name 
                self.out + '( max_length * ' + max + ')'
                self.out + ')'
        
    def output_runloop_function_def_start(self):
        self.out.lf() + 'subroutine run_loop'
        self.out.indent()
        self.out.lf() + 'use iso_c_binding'
        self.out.lf() + 'use forsockets'
        self.output_modules()
        self.out.lf() + "implicit none"
        self.out.lf()
        self.output_mpi_include()
        self.out.lf().lf()
        self.out.lf().lf() + 'integer (c_int32_t) internal__redirect_outputs'
        self.out.lf().lf()
        self.out.n() + 'integer (c_int32_t) :: max_length = 255, MAX_STRING_LENGTH = ' + self.MAX_STRING_LEN
        self.out.n() + 'logical :: must_run_loop, error'
        self.out.n() + 'integer i, length'
        
        self.out.lf() + 'character (c_char), allocatable, target :: characters_in(:) * ' + self.MAX_STRING_LEN
        self.out.lf() + 'character (c_char), allocatable, target :: characters_out(:) * ' + self.MAX_STRING_LEN
        
        #self.out.n() + 'character (len=100000) :: output_characters'
        #self.out.lf().lf() + 'character (c_char), target :: input_characters(100)'
        
        self.out.lf().lf() + 'integer (c_int32_t), target :: header_in(HEADER_SIZE)'
        self.out.n() + 'integer (c_int32_t), target :: header_out(HEADER_SIZE)'
        self.out.lf()
        
        self.output_declarations_for_the_functions()
        self.out.lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + spec.type
            self.out + ', allocatable, target :: '
            self.out + spec.input_var_name
            self.out + '(:)'
            
            self.out.lf() + spec.type
            self.out + ', allocatable, target :: '
            self.out + spec.output_var_name
            self.out + '(:)'
                        
#            self.out.lf() + 'integer ::' + ' '
#            self.out + dtype_spec.counter_name + '_out'
#            self.out + ', ' + dtype_spec.counter_name + '_in'
            self.out.lf()
            
        
        self.out.lf()
        self.output_allocate_arrays()
#        self.out.lf()
#        self.out.lf().lf() + 'call MPI_COMM_GET_PARENT(parent, ioerror)'
#        self.out.lf()      + 'call MPI_COMM_RANK(parent, rank, ioerror)'
        self.out.lf().lf() + 'must_run_loop = .true.'
        self.out.lf().lf() + 'do while (must_run_loop)'
        self.out.indent()

        self.out.lf().lf() + 'call receive_integers(c_loc(header_in), HEADER_SIZE)'
        self.out.lf().lf() + 'length = header_in(HEADER_CALL_COUNT)'
        
        self.out.lf().lf() + 'if (length .gt. max_length) then'
        self.out.indent()
        self.out.lf() + 'max_length = length + 255;'
        self.output_deallocate_statements()
        self.output_allocate_arrays()
        self.out.dedent()
        self.out.lf() + 'end if'
        self.out.lf()
    
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'if (header_in(' + spec.counter_name + ')'
            self.out + ' .gt. 0) then'
    
            self.out.indent().lf()
            
            if dtype == 'string':
                self.out.lf() + string_receive_code
            else:
                self.out + 'call receive_'
                self.out + spec.mpi_type + 's(c_loc('
                self.out + spec.input_var_name + '), header_in('
                self.out + spec.counter_name + '))'
        
            self.out.dedent().lf()
            self.out + 'end if'
            self.out.lf()
            
            
        self.out.lf().lf() + 'header_out = 0'
        self.out.lf() + 'header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)'
        self.out.lf() + 'header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)'
        self.out.lf() + 'header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)'
        
        self.out.lf().lf() + 'error = .false.'
        
         
    def output_switch_start(self):
        self.out.lf().lf() + 'select case (header_in(HEADER_FUNCTION_ID))'
        self.out.indent()
        self.out.lf() + 'case(0)'
        self.out.indent().lf()+'must_run_loop = .false.'
        self.out.dedent()
        
    def output_switch_end(self):
        self.out.lf() + 'case default'
        self.out.indent().lf() + 'error = .true.'
        self.out.dedent()
        self.out.dedent().lf() + 'end select'
        
    def output_runloop_function_def_end(self):
        
#        self.out.lf().lf() + 'if (rank .eq. 0 ) then'
#        self.out.indent().lf()
#        
#        self.out.lf().lf() + 'header(1) = tag_out'
#        self.out.lf() + 'header(2) = len_out'
#        for i, dtype in enumerate(dtypes):
#            spec = self.dtype_to_spec[dtype]
#            self.out.lf()  + 'header(' + (i+3) + ')'
#            self.out + ' = ' + spec.counter_name + '_out'
#        
#        self.out.lf().lf() + 'call MPI_SEND(header, '
#        self.out + self.length_of_the_header
#        self.out + ', MPI_INTEGER,'
#        self.out + ' 0, 999, &'
#        self.out.indent().lf() + 'parent, ioerror);'
#        self.out.dedent().lf()
        
        self.out.lf().lf() + 'if (header_in(HEADER_STRING_COUNT) .gt. 0) then'
        self.out.indent()
        self.out.lf().lf() + 'deallocate(characters_in)'
        self.out.dedent()
        self.out.lf().lf() + 'end if'
   
        self.out.lf().lf() + "!print*, 'sending header', header_out"
        self.out.lf().lf() + '!call flush()'
        self.out.lf().lf() + 'call send_integers(c_loc(header_out), HEADER_SIZE)'
        self.out.lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            
            self.out.lf() + 'if (header_out(' + spec.counter_name + ')'
            self.out + ' .gt. 0) then'
            
            self.out.indent()
            
            if dtype == 'string':
                self.out.lf() + string_send_code
                self.out.lf().lf() + 'deallocate(characters_out)'
            else:
                self.out.lf() + 'call send_'
                self.out + spec.mpi_type + 's(c_loc('
                self.out + spec.output_var_name + '), header_out('
                self.out + spec.counter_name + '))'
                
            self.out.dedent().lf() +'end if'
            self.out.lf()
        
        self.out.dedent()
        self.out.lf() + 'end do'
        
        
        self.out.lf()
        self.output_deallocate_statements()
            
#        self.out.lf() + 'call MPI_COMM_DISCONNECT(parent, ioerror)'
        self.out.lf() + 'return'
        
        self.out.dedent().lf() + 'end subroutine'
    
    def output_deallocate_statements(self):
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'deallocate(' + dtype_spec.input_var_name  + ')'
            
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            if max > 0:
                self.out.lf() + 'deallocate(' + dtype_spec.output_var_name  + ')'
    
    def output_forsockets_module(self):
        self.out + forsockets_module_code
        
    def output_redirect_output(self):
        self.out.lf().lf() + redirect_outputs_function_code
        
    def output_main(self):
        self.out.lf().lf() + main_program_code
        
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
        
    
        
        

        
       
    
        
        
        
