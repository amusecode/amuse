from amuse.support.core import late
from amuse.legacy.support.core import RemoteFunction
from amuse.legacy.support.create_code import MakeCodeString
from amuse.legacy.support.create_code import MakeCodeStringOfAClassWithLegacyFunctions
from amuse.legacy.support.create_code import DTypeSpec
        
        
class MakeAFortranStringOfALegacyFunctionSpecification(MakeCodeString):
    
    @late
    def dtype_to_spec(self):
        return {
            'i' : DTypeSpec('integers_in','integers_out', 
                            'number_of_integers_out', 'integer'),
            'd' : DTypeSpec('doubles_in', 'doubles_out',
                            'number_of_doubles_out', 'real*8'),
            'f' : DTypeSpec('floats_in', 'floats_out',
                            'number_of_floats_out', 'real*4')}
   
    def start(self):
        
        if not self.specification.result_type is None:
            raise Exception("cannot handle result type for "\
            "fortran subroutines! (no support for fortran "\
            "functions yet)")
            
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
        
        for name, dtype, direction in self.specification.parameters:
            spec = self.dtype_to_spec[dtype]
            
            if first:
                first = False
                self.out + ' &'
            else:
                self.out + ' ,&'
                
            if direction == RemoteFunction.IN:
                self.out.n() + spec.input_var_name 
                self.out + '(' + (spec.number_of_inputs + 1) + ')'
                spec.number_of_inputs += 1
            if direction == RemoteFunction.INOUT:
                self.out.n() + spec.input_var_name 
                self.out + '(' + (spec.number_of_inputs + 1) + ')'
                spec.number_of_inputs += 1
            elif direction == RemoteFunction.OUT:
                self.out.n() + spec.output_var_name
                self.out + '(' + (spec.number_of_outputs + 1) + ')'
                spec.number_of_outputs += 1
    
        self.out.dedent()
        
    def output_lines_with_inout_variables(self):
        dtype_to_incount = {}
        
        for name, dtype, direction in self.specification.parameters:
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_incount.get(dtype, 0)
        
            if direction == RemoteFunction.IN:
                dtype_to_incount[dtype] = count + 1
            if direction == RemoteFunction.INOUT:
                self.out.n() + spec.output_var_name 
                self.out + '(' + (spec.number_of_outputs + 1)  + ')' 
                self.out + ' = ' 
                self.out + spec.input_var_name + '(' + (count + 1) + ')'
                spec.number_of_outputs += 1
                dtype_to_incount[dtype] = count + 1
    
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT \
                or direction == RemoteFunction.INOUT:
                count = dtype_to_count.get(dtype, 0)
                dtype_to_count[dtype] = count + 1
                
        if not self.specification.result_type is None:
            count = dtype_to_count.get(self.specification.result_type, 0)
            dtype_to_count[self.specification.result_type] = count + 1
            
        for dtype in dtype_to_count:       
            spec = self.dtype_to_spec[dtype]
            count = dtype_to_count[dtype]
            self.out.n() 
            self.out + spec.output_counter_name + ' = ' + count 
            pass
            
    def output_function_end(self):
        self.out + ' &'
        self.out.n() + ')'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.output_var_name
            self.out + '(' + (spec.number_of_outputs + 1) + ')' + ' = '
            spec.number_of_outputs += 1
        self.out + 'CALL ' +  self.specification.name + '('
        
    def output_casestmt_start(self):
        self.out + 'CASE(' + self.specification.id + ')'
        
    def output_casestmt_end(self):
        self.out.n() 
        
        
        

class MakeAFortranStringOfAClassWithLegacyFunctions \
    (MakeCodeStringOfAClassWithLegacyFunctions):

    @late
    def dtype_to_spec(self):
        return {
            'i' : DTypeSpec('integers_in','integers_out',
                            'number_of_integers', 'integer'),
            'd' : DTypeSpec('doubles_in', 'doubles_out',
                            'number_of_doubles', 'real*8'),
            'f' : DTypeSpec('floats_in', 'floats_out',
                            'number_of_floats_out', 'real*4')}
   
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
        self.out.lf().lf() + 'integer header(4)'
        self.out.lf().lf() + 'integer :: tag_in, tag_out'
        
        maximum_number_of_inputvariables_of_a_type = 255
        for dtype_spec in self.dtype_to_spec.values():
            
            self.out.lf() + dtype_spec.type + ' ' 
            self.out + dtype_spec.input_var_name 
            self.out + '(' + maximum_number_of_inputvariables_of_a_type + ')'
            
            self.out.lf() + dtype_spec.type 
            self.out + ' ' + dtype_spec.output_var_name 
            self.out + '(' + maximum_number_of_inputvariables_of_a_type + ')'
            
            self.out.lf() + 'integer ::' + ' ' 
            self.out + dtype_spec.output_counter_name + '_out' 
            self.out + ', ' + dtype_spec.output_counter_name + '_in'
            
        self.out.lf().lf() + 'call MPI_COMM_GET_PARENT(parent, ioerror)'
        self.out.lf()      + 'call MPI_COMM_RANK(parent, rank, mpierror)'
        self.out.lf().lf() + 'must_run_loop = 1'
        self.out.lf().lf() + 'do while (must_run_loop .eq. 1)'
        self.out.indent()
       
       
        self.out.lf() + 'call MPI_RECV(header, 4, MPI_INTEGER, 0,'
        self.out + ' rank, parent,&'
        self.out.indent().lf() + 'mpiStatus, ioerror)'
        self.out.dedent()
        self.out.lf().lf() + 'tag_in = header(1)'
        
        spec = [
          ('number_of_doubles_in', 'doubles_in', 'MPI_DOUBLE_PRECISION')
          ,('number_of_integers_in', 'integers_in', 'MPI_INTEGER')
          ,('number_of_floats_in', 'floats_in', 'MPI_SINGLE_PRECISION')]
        
        for i, (number_parameter, input_parameter_name, mpi_type)\
            in enumerate(spec):
            self.out.lf() + number_parameter + ' =  ' 
            self.out + 'header(' + (i+2) + ')'
        
        self.out.lf().lf() + 'tag_out = tag_in'
        
        
        self.out.lf()      + 'number_of_doubles_out = 0'
        self.out.lf()      + 'number_of_integers_out = 0'
        self.out.lf()      + 'number_of_floats_out = 0'
        self.out.lf()
        
       
        for number_parameter, input_parameter_name, mpi_type in spec:
               self.out.lf() + 'if (' + number_parameter + ' .gt. 0) then'
               
               self.out.indent().lf() + 'call MPI_RECV('
               self.out + input_parameter_name + ', ' 
               self.out + number_parameter + ', &'
               
               self.out.indent().n() + mpi_type+ ', 0, rank, parent,&'
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
        self.out.lf()      + 'header(2) = number_of_doubles_out'
        self.out.lf()      + 'header(3) = number_of_integers_out'
        
        self.out.lf().lf() + 'call MPI_SEND(header, 4, MPI_INTEGER,'
        self.out + ' 0, 999, &'
        self.out.indent().lf() + 'parent, mpierror);'
        self.out.dedent().lf()
        
        spec = [
          ('number_of_doubles_out', 'doubles_out', 'MPI_DOUBLE_PRECISION')
          ,('number_of_integers_out', 'integers_out', 'MPI_INTEGER')
          ,('number_of_floats_out', 'floats_out', 'MPI_SINGLE_PRECISION')]
        for number_parameter, parameter_name, mpi_type in spec:
               self.out.lf() + 'if (' + number_parameter + ' .gt. 0) then'
               self.out.indent().lf() 
               self.out + 'call MPI_SEND(' 
               self.out + parameter_name + ', ' + number_parameter + ', &'
               self.out.indent().lf() + mpi_type + ', 0, 999, &'
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
        
       
    
        
        
        
