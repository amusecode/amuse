from amuse.support.core import late
from amuse.legacy.support.core import RemoteFunction
from amuse.legacy.support.create_code import MakeCodeString
from amuse.legacy.support.create_code import MakeCodeStringOfAClassWithLegacyFunctions
from amuse.legacy.support.create_code import DTypeSpec, dtypes
 
dtype_to_spec = {
    'i' : DTypeSpec('integers_in','integers_out', 
                    'number_of_integers', 'integer', 'MPI_INTEGER'),
    'd' : DTypeSpec('doubles_in', 'doubles_out',
                    'number_of_doubles', 'real*8', 'MPI_DOUBLE_PRECISION'),
    'f' : DTypeSpec('floats_in', 'floats_out',
                    'number_of_floats', 'real*4', 'MPI_SINGLE_PRECISION')} 
        
class MakeAFortranStringOfALegacyFunctionSpecification(MakeCodeString):
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
   
    def start(self):
                   
        self.specification.prepare_output_parameters()
         
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
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.dtype]
            
            if first:
                first = False
                self.out + ' &'
            else:
                self.out + ' ,&'
                
            if parameter.direction == RemoteFunction.IN:
                self.out.n() + spec.input_var_name 
                self.out + '(' + (parameter.input_index + 1) + ')'
            if parameter.direction == RemoteFunction.INOUT:
                self.out.n() + spec.input_var_name 
                self.out + '(' + (parameter.input_index + 1) + ')'
            elif parameter.direction == RemoteFunction.OUT:
                self.out.n() + spec.output_var_name
                self.out + '(' + (parameter.output_index + 1) + ')'
                
        self.out.dedent()
        
    def output_lines_with_inout_variables(self):
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.dtype]
            
            if parameter.direction == RemoteFunction.INOUT:
                self.out.n() + spec.output_var_name 
                self.out + '(' + (parameter.output_index + 1)  + ')' 
                self.out + ' = ' 
                self.out + spec.input_var_name + '(' + (parameter.input_index  + 1) + ')'
    
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for parameter in self.specification.output_parameters:
            count = dtype_to_count.get(parameter.dtype, 0)
            dtype_to_count[parameter.dtype] = count + 1
                
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
            self.out + spec.output_var_name
            self.out + '(' + 1 + ')' + ' = '
        else:    
            self.out + 'CALL ' 
        self.out +  self.specification.name + '('
        
    def output_casestmt_start(self):
        self.out + 'CASE(' + self.specification.id + ')'
        
    def output_casestmt_end(self):
        self.out.n() 
        
        
class MakeAFortranStringOfALegacyGlobalSpecification(MakeCodeString):
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
            
    def start(self):
        self.output_casestmt_start()
        self.out.indent()
        
        spec = self.dtype_to_spec[self.legacy_global.dtype]
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

class MakeAFortranStringOfAClassWithLegacyFunctions \
    (MakeCodeStringOfAClassWithLegacyFunctions):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec 
   
    @late
    def number_of_types(self):
        return len(self.dtype_to_spec)
        
    @late
    def length_of_the_header(self):
        return 2 + self.number_of_types
        
        
    def make_legacy_function(self):
        return MakeAFortranStringOfALegacyFunctionSpecification()
   
    def start(self):
        self.output_runloop_function_def_start()
        self.output_switch_start()
        self.output_legacy_functions()
        self.output_switch_end()
        self.output_runloop_function_def_end()
        self.output_main()
        self._result = self.out.string

    def output_mpi_include(self):
        self.out.n() + "INCLUDE 'mpif.h'"
        
    def output_local_includes(self):
        self.out.n()
        for x in ['muse_dynamics.h', 'parameters.h', 'local.h']:
            self.out.n() + '#include "' + x + '"'
        
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'SUBROUTINE run_loop'
        self.out.indent()
        self.output_mpi_include()
        self.out.n() + 'integer :: rank, parent, ioerror'
        self.out.n() + 'integer :: must_run_loop'
        self.out.n() + 'integer mpiStatus(MPI_STATUS_SIZE,4)'
        self.out.lf().lf() + 'integer header(' 
        self.out + self.length_of_the_header + ')'
        self.out.lf().lf() + 'integer :: tag_in, tag_out'
        self.out.lf().lf() + 'integer :: len_in, len_out'
        
        maximum_number_of_inputvariables_of_a_type = 255
        for dtype_spec in self.dtype_to_spec.values():
            
            self.out.lf() + dtype_spec.type + ' ' 
            self.out + dtype_spec.input_var_name 
            self.out + '(' + maximum_number_of_inputvariables_of_a_type + ')'
            
            self.out.lf() + dtype_spec.type 
            self.out + ' ' + dtype_spec.output_var_name 
            self.out + '(' + maximum_number_of_inputvariables_of_a_type + ')'
            
            self.out.lf() + 'integer ::' + ' ' 
            self.out + dtype_spec.counter_name + '_out' 
            self.out + ', ' + dtype_spec.counter_name + '_in'
            
        self.out.lf().lf() + 'call MPI_COMM_GET_PARENT(parent, ioerror)'
        self.out.lf()      + 'call MPI_COMM_RANK(parent, rank, mpierror)'
        self.out.lf().lf() + 'must_run_loop = 1'
        self.out.lf().lf() + 'do while (must_run_loop .eq. 1)'
        self.out.indent()
       
       
        self.out.lf() + 'call MPI_RECV(header, '
        self.out + self.length_of_the_header
        self.out + ', MPI_INTEGER, 0,'
        self.out + ' 0, parent,&'
        self.out.indent().lf() + 'mpiStatus, ioerror)'
        self.out.dedent()
        self.out.lf().lf() + 'tag_in = header(1)'
        self.out.lf().lf() + 'len_in = header(2)'
        
        
               
               
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + spec.counter_name + '_in' + ' =  ' 
            self.out + 'header(' + (i+3) + ')'
        
        self.out.lf().lf() + 'tag_out = tag_in'
        self.out.lf() + 'len_out = len_in'
        
        
        self.out.lf()      + 'number_of_doubles_out = 0'
        self.out.lf()      + 'number_of_integers_out = 0'
        self.out.lf()      + 'number_of_floats_out = 0'
        self.out.lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'if (' + spec.counter_name + '_in' 
            self.out + ' .gt. 0) then'

            self.out.indent().lf() + 'call MPI_RECV('
            self.out + spec.input_var_name + ', ' 
            self.out + spec.counter_name + '_in' + ', &'

            self.out.indent().n() + spec.mpi_type
            self.out + ', 0, 0, parent,&'
            self.out.n() + 'mpiStatus, ioError);'

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
        
        self.out.lf().lf() + 'header(1) = tag_out'
        self.out.lf() + 'header(2) = len_out'
        for i, dtype in enumerate(['d','i','f']):
            spec = self.dtype_to_spec[dtype]
            self.out.lf()  + 'header(' + (i+3) + ')'
            self.out + ' = ' + spec.counter_name + '_out'
        
        self.out.lf().lf() + 'call MPI_SEND(header, '
        self.out + self.length_of_the_header 
        self.out + ', MPI_INTEGER,'
        self.out + ' 0, 999, &'
        self.out.indent().lf() + 'parent, mpierror);'
        self.out.dedent().lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'if (' + spec.counter_name + '_out'
            self.out + ' .gt. 0) then'
            self.out.indent().lf() 
            self.out + 'call MPI_SEND(' 
            self.out + spec.output_var_name + ', ' +  spec.counter_name + '_out'
            self.out + ', &'
            self.out.indent().lf() + spec.mpi_type + ', 0, 999, &'
            self.out.lf() + 'parent, mpierror);'
            self.out.dedent().dedent().lf() +'end if'
        self.out.dedent()
        self.out.lf() + 'end do'
        self.out.lf() + 'return'
        
        self.out.dedent().lf() + 'end subroutine'
        
    def output_main(self):
        self.out.lf().lf() + 'program muse_worker'
        self.out.indent()
        self.output_mpi_include()
        self.out.lf() + 'call MPI_INIT(mpierror)'
        self.out.lf().lf() + 'call run_loop()'
        self.out.lf().lf() + 'call MPI_FINALIZE(mpierror)'
        self.out.dedent().lf()+'end program muse_worker'
        
       
    
        
        
        
