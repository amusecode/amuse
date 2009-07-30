from amuse.legacy.support.core import *
        
        
class MakeAFortranStringOfALegacyFunctionSpecification(object):
    def __init__(self):
        pass
    
    @late  
    def result(self):
        self.start()
        return self._result
    
    @late
    def dtype_to_spec(self):
        class DTypeSpec(object):
            def __init__(self, input_var_name, output_var_name, output_counter_name):
                self.number_of_inputs = 0
                self.number_of_outputs = 0
                self.input_var_name = input_var_name
                self.output_var_name = output_var_name
                self.output_counter_name = output_counter_name
             
        return {
            'i' : DTypeSpec('integers_in', 'integers_out', 'number_of_integers_out'),
            'd' : DTypeSpec('doubles_in', 'doubles_out', 'number_of_doubles_out')}
    @late
    def out(self):
        return print_out()
                
    def start(self):
        
        if not self.specification.result_type is None:
            raise Exception("cannot handle result type for fortran subroutines!  (no support for fortran functions yet)")
            
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
                self.out.n() + spec.input_var_name + '(' + (spec.number_of_inputs + 1) + ')'
                spec.number_of_inputs += 1
            if direction == RemoteFunction.INOUT:
                self.out.n() + spec.input_var_name + '(' + (spec.number_of_inputs + 1) + ')'
                spec.number_of_inputs += 1
            elif direction == RemoteFunction.OUT:
                self.out.n() + spec.output_var_name + '(' + (spec.number_of_outputs + 1) + ')'
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
                self.out.n() + spec.output_var_name + '(' + (spec.number_of_outputs + 1)  + ')' + ' = ' + spec.input_var_name + '(' + (count + 1) + ')'
                spec.number_of_outputs += 1
                dtype_to_incount[dtype] = count + 1
    
    def output_lines_with_number_of_outputs(self):
        dtype_to_count = {}
        
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT or direction == RemoteFunction.INOUT:
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
            self.out + spec.c_output_var_name+ '(' + (spec.number_of_outputs + 1) + ')' + ' = '
            spec.number_of_outputs += 1
        self.out + 'CALL ' +  self.specification.name + '('
        
    def output_casestmt_start(self):
        self.out + 'CASE(' + self.specification.id + ')'
        
    def output_casestmt_end(self):
        self.out.n() 
        
        
        
HEADER_CLASS_STRING = """class message_header {

public:
	int tag;
	int number_of_doubles;
	int number_of_ints;
	message_header(): tag(0), number_of_doubles(0), number_of_ints(0) {}

	void send(MPI::Intercomm & intercomm, int rank){
		int header[3];
		header[0] = tag;
		header[1] = number_of_doubles;
		header[2] = number_of_ints;		
		intercomm.Send(header, 3, MPI_INT, 0, 999);
	}

	void recv(MPI::Intercomm & intercom, int rank) {
		int header[6];

		intercom.Recv(header, 3, MPI_INT, 0, rank);
		tag = header[0];
		number_of_doubles = header[1];
		number_of_ints = header[2];
	}
};"""
        

class MakeAFortranStringOfAClassWithLegacyFunctions(object):
    def __init__(self):
        pass
    
    @late  
    def result(self):
        self.start()
        return self._result
    
    @late
    def dtype_to_spec(self):
        class DTypeSpec(object):
            def __init__(self, c_input_var_name, c_output_var_name, c_output_counter_name, c_type):
                self.number_of_inputs = 0
                self.number_of_outputs = 0
                self.c_input_var_name = c_input_var_name
                self.c_output_var_name = c_output_var_name
                self.c_output_counter_name = c_output_counter_name
                self.c_type = c_type
             
        return {
            'i' : DTypeSpec('ints_in','ints_out', 'number_of_ints', 'int'),
            'd' : DTypeSpec('doubles_in', 'doubles_out','number_of_doubles', 'double')}
    @late
    def out(self):
        return print_out()
        
    def start(self):
        self.output_mpi_include()
        self.output_local_includes()
        self.out.lf().lf()
        self.output_header_class_definition()
        self.output_runloop_function_def_start()
        self.output_switch_start()
        
        attribute_names = dir(self.class_with_legacy_functions)
        legacy_functions = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(self.class_with_legacy_functions, x)
            if isinstance(value, legacy_function):
                legacy_functions.append(value)
                
        for x in legacy_functions:
            if x.specification.id == 0:
                continue
            self.out.lf()
            uc = MakeACStringOfALegacyFunctionSpecification()
            uc.specification = x.specification
            uc.out = self.out
            uc.start()
                 
        self.output_switch_end()
        self.output_runloop_function_def_end()
        self.output_main()
        self._result = self.out.string
        
    def output_mpi_include(self):
        self.out.n() + '#include <mpi.h>'
    def output_local_includes(self):
        self.out.n()
        for x in ['muse_dynamics.h', 'parameters.h', 'local.h']:
            self.out.n() + '#include "' + x + '"'
    def output_header_class_definition(self):
        self.out + HEADER_CLASS_STRING
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'void run_loop() {'
        self.out.indent()
        self.out.n() + 'int rank = MPI::COMM_WORLD.Get_rank();'
        self.out.lf().lf() + 'MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();'
        self.out.lf().lf() + 'bool must_run_loop = true;'
        self.out.lf().lf() + 'while(must_run_loop) {'
        self.out.indent()
        maximum_number_of_inputvariables_of_a_type = 255
        for dtype_spec in self.dtype_to_spec.values():
            self.out.lf() + dtype_spec.c_type + ' ' + dtype_spec.c_input_var_name + '[' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
            self.out.lf() + dtype_spec.c_type + ' ' + dtype_spec.c_output_var_name + '[' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
        self.out.lf()
        self.out.lf() + 'message_header request_header;'
        self.out.lf() + 'message_header reply_header;'
        self.out.lf()
        self.out.lf() + 'request_header.recv(parent,rank);'
        spec = [('number_of_doubles', 'doubles_in', 'MPI_DOUBLE'),('number_of_ints', 'ints_in', 'MPI_INT')]
        for number_parameter, input_parameter_name, mpi_type in spec:
               self.out.lf() + 'if(request_header.' + number_parameter + ' > 0) {'
               self.out.indent().lf() + 'parent.Recv(' + input_parameter_name + ', ' + 'request_header.' + number_parameter + ', ' + mpi_type+ ', 0, rank);'
               self.out.dedent().lf() +'}'
        self.out.lf().lf() + 'reply_header.tag = request_header.tag;'
    def output_switch_start(self):
        self.out.lf().lf() + 'switch(request_header.tag) {'
        self.out.indent()
        self.out.lf() + 'case 0:'
        self.out.indent().lf()+'must_run_loop = false;'
        self.out.lf()+'break;'
        self.out.dedent()
    def output_switch_end(self):
        self.out.lf() + 'default:'
        self.out.indent().lf() + 'reply_header.tag = -1;'
        self.out.dedent()
        self.out.dedent().lf() + '}'
        
        
    
            
            
        
        
    def output_runloop_function_def_end(self):
        self.out.lf().lf() + 'reply_header.send(parent, rank);'
        spec = [('number_of_doubles', 'doubles_out', 'MPI_DOUBLE'),('number_of_ints', 'ints_out', 'MPI_INT')]
        for number_parameter, parameter_name, mpi_type in spec:
               self.out.lf() + 'if(reply_header.' + number_parameter + ' > 0) {'
               self.out.indent().lf() + 'parent.Send(' + parameter_name + ', ' + 'reply_header.' + number_parameter + ', ' + mpi_type+ ', 0, 999);'
               self.out.dedent().lf() +'}'
        self.out.dedent()
        self.out.lf() + '}'
        self.out.dedent()
        self.out.lf() + '}'
    def output_main(self):
        self.out.lf().lf() + 'int main(int argc, char *argv[])'
        self.out.lf() + '{'
        self.out.indent().lf() + 'MPI::Init(argc, argv);'
        self.out.lf().lf() + 'run_loop();'
        self.out.lf().lf() + 'MPI_Finalize();'
        self.out.lf() + 'return 0;'
        self.out.dedent().lf()+'}'
    
        
        
        
