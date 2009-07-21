"""
"""
from amuse.support.core import late
from amuse.support.core import print_out

from zlib import crc32

import numpy

class LegacyCall(object):
    """A legacy_call implements the runtime call to the remote process.
    """
    def __init__(self, interface, owner, specification):
        self.interface = interface
        self.owner = owner
        self.specification = specification
    
    def __call__(self, *arguments_list, **keyword_arguments):
        dtype_to_values_and_keyword = {
            'd' : [[],'doubles_in',0],
            'i' : [[],'ints_in',0]
        }
        
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.IN or direction == RemoteFunction.INOUT:
                values = dtype_to_values_and_keyword[dtype]
                values[0].append(None)
            
        names_in_argument_list = set([])
        for index, argument in enumerate(arguments_list):
            name, dtype, direction = self.specification.parameters[index]
            names_in_argument_list.add(name)
            values = dtype_to_values_and_keyword[dtype]
            values[0][values[2]] = argument
            values[2] += 1
        
        for parameters in self.specification.dtype_to_input_parameters.values():
            for index, parameter in enumerate(parameters):
                name, dtype, direction = parameter
                if name in keyword_arguments:
                    values = dtype_to_values_and_keyword[dtype]
                    values[0][index] = keyword_arguments[name]
               
        call_keyword_arguments = {}
        for values, keyword, count in dtype_to_values_and_keyword.values():
            call_keyword_arguments[keyword] = values
            
        return self.do_call(self.specification.id, **call_keyword_arguments)
        
    def do_call(self, id, **keyword_arguments):
        self.interface.channel.send_message(id , **keyword_arguments)
        (doubles, ints) = self.interface.channel.recv_message(id)
        
        dtype_to_values_and_keyword = {
            'd' : [[],'doubles_in',0],
            'i' : [[],'ints_in',0]
        }
        
        number_of_outputs = 0
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT or direction == RemoteFunction.INOUT:
                values = dtype_to_values_and_keyword[dtype]
                values[0].append(None)
                number_of_outputs += 1
       
        
        if number_of_outputs == 0:
            if self.specification.result_type == 'i':
                return ints[0]       
            if self.specification.result_type == 'd':
                return doubles[0]       
        
        return (doubles, ints)
       

class legacy_function(object):
    """The meta information for a function call to a code
    """
    def __init__(self, specification_function):
        self.specification_function = specification_function
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return LegacyCall(instance, owner, self.specification)
        
    def __set__(self, instance, value):
        return
        
        
    
    def to_c_string(self):
        uc = MakeACStringOfALegacyFunctionSpecification()
        uc.specification = self.specification
        return uc.result 
        
    @late
    def specification(self):
        result = self.specification_function()
        if result.name is None:
            result.name = self.specification_function.__name__
        if result.id is None:
            result.id = crc32(result.name)
        return result
        
        
class legacy_global(object):
    """The meta information for globals of the code
    """
    def __init__(self, name , id = None, dtype = 'i'):
        self.name = name
        self.id = id
        self.dtype = dtype
        
        if self.id is None:
            self.id = crc32(self.name)
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return LegacyCall(instance, owner, self.get_specification)()
        
    def __set__(self, instance, value):
        return LegacyCall(instance, None, self.set_specification)(value)
        
        
    
    def to_c_string(self):
        uc = MakeACStringOfALegacyGlobalSpecification()
        uc.specification = self
        return uc.result 
        
    @late
    def set_specification(self):
        result = RemoteFunction()
        result.id = self.id
        result.name = self.name
        result.addParameter('value', dtype=self.dtype, direction=RemoteFunction.IN)
        return result
    @late
    def get_specification(self):
        result = RemoteFunction()
        result.id = self.id
        result.name = self.name
        result.result_type = self.dtype
        return result
        
class RemoteFunction(object):
    IN = object()
    OUT = object()
    INOUT = object()
    
    def __init__(self):
        self.parameters = []
        self.name = None
        self.id = None
        self.result_type = None
        
    def addParameter(self, name, dtype = 'i', direction = IN):
        self.parameters.append((name, dtype, direction))
    
    @property
    def dtype_to_input_parameters(self):
        result = {}
        for name, dtype, direction in self.parameters:
            parameters = result.get(dtype, [])
            parameters.append((name, dtype, direction))
            result[dtype] = parameters
        return result
                
    
    
        
                
        
        
class MakeACStringOfALegacyFunctionSpecification(object):
    def __init__(self):
        pass
    
    @late  
    def result(self):
        self.start()
        return self._result
    
    @late
    def dtype_to_spec(self):
        class DTypeSpec(object):
            def __init__(self, c_input_var_name, c_output_var_name, c_output_counter_name):
                self.number_of_inputs = 0
                self.number_of_outputs = 0
                self.c_input_var_name = c_input_var_name
                self.c_output_var_name = c_output_var_name
                self.c_output_counter_name = c_output_counter_name
             
        return {
            'i' : DTypeSpec('ints_in','ints_out', 'number_of_ints'),
            'd' : DTypeSpec('doubles_in', 'doubles_out','number_of_doubles')}
    @late
    def out(self):
        return print_out()
        
        
    def start(self):
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
            else:
                self.out + ' ,'
                
            if direction == RemoteFunction.IN:
                self.out.n() + spec.c_input_var_name + '[' + spec.number_of_inputs + ']'
                spec.number_of_inputs += 1
            if direction == RemoteFunction.INOUT:
                self.out.n() + '&' + spec.c_input_var_name + '[' + spec.number_of_inputs + ']'
                spec.number_of_inputs += 1
            elif direction == RemoteFunction.OUT:
                self.out.n() + '&' + spec.c_output_var_name + '[' + spec.number_of_outputs + ']'
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
                self.out.n() + spec.c_output_var_name + '[' + spec.number_of_outputs + ']' + ' = ' + spec.c_input_var_name + '[' + count + ']'+';'
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
            self.out + 'reply.' + spec.c_output_counter_name + ' = ' + count + ';'
            pass
    def output_function_end(self):
        self.out.n() + ')' + ';'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.c_output_var_name+ '[' + spec.number_of_outputs + ']' + ' = '
            spec.number_of_outputs += 1
        self.out + self.specification.name + '('
    def output_casestmt_start(self):
        self.out + 'case ' + self.specification.id + ':'
    def output_casestmt_end(self):
        self.out.n() + 'break;'
        
        
        
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
        

class MakeACStringOfAClassWithLegacyFunctions(object):
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
    
        
        
        
class MpiChannel(object):
    from mpi4py import MPI
    
    def __init__(self, intercomm):
        self.intercomm = intercomm
        
    def send_message(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_in=[], ints_in=[]):
        header = numpy.array([tag,len(doubles_in), len(ints_in)], dtype='i')
        self.intercomm.Send([header, self.MPI.INT], dest=0, tag=0)
        if doubles_in:
            doubles = numpy.array(doubles_in, dtype='d')
            self.intercomm.Send([doubles, self.MPI.DOUBLE], dest=0, tag=0)
        if ints_in:
            ints = numpy.array(ints_in, dtype='i')
            self.intercomm.Send([ints, self.MPI.INT], dest=0, tag=0)
            
    def recv_message(self, tag):
        header = numpy.empty(3,  dtype='i')
        self.intercomm.Recv([header, self.MPI.INT], source=0, tag=999)
        n_doubles = header[1]
        n_ints = header[2]
        if n_doubles > 0:
            doubles_result = numpy.empty(n_doubles,  dtype='d')
            self.intercomm.Recv([doubles_result, self.MPI.DOUBLE], source=0, tag=999)
        else:
            doubles_result = []
        if n_ints > 0:
            ints_result = numpy.empty(n_ints,  dtype='i')
            self.intercomm.Recv([ints_result, self.MPI.INT], source=0, tag=999)
        else:
            ints_result = []
        if header[0] < 0:
            raise Exception("Not a valid message!")
        return (doubles_result, ints_result)
    


        
            