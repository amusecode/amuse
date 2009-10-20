from amuse.support.core import late
from amuse.legacy.support.core import RemoteFunction
from amuse.legacy.support.create_code import MakeCodeString
from amuse.legacy.support.create_code import MakeCodeStringOfAClassWithLegacyFunctions
from amuse.legacy.support.create_code import DTypeSpec, dtypes, DTypeToSpecDictionary
      
import numpy

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('ints_in','ints_out',
                    'number_of_ints', 'int', 'MPI_INT'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'number_of_doubles', 'double', 'MPI_DOUBLE'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'number_of_floats', 'float', 'MPI_FLOAT'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'number_of_strings', 'int', 'MPI_INTEGER'),
})


class MakeCCodeString(MakeCodeString):
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
       
         
class MakeACStringOfALegacyFunctionSpecification(MakeCCodeString):

        
    def start(self):
        
        self.specification.prepare_output_parameters()
        self.output_casestmt_start()
        self.out.indent()
        
        if self.specification.can_handle_array:
            self.out.lf() + 'for (int i = 0 ; i < request_header.len; i++){'
            self.out.indent()
        
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        
        if self.specification.can_handle_array:
            self.out.dedent()
            self.out.lf() + '}'
        
        self.output_lines_with_number_of_outputs()
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
    
    def index_string(self, index):
        if self.specification.can_handle_array:
            if index == 0:
                return 'i'
            else:
                return '( %d * request_header.len) + i' % index
        else:
            return index
            
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ' ,'
                
            if parameter.direction == RemoteFunction.IN:
                if parameter.datatype == 'string':
                    self.out.n() + 'characters + ' 
                    self.out + '( ' + self.index_string(parameter.input_index) + '- 1 < 0 ? 0 :' 
                    self.out + spec.input_var_name
                    self.out + '[' + self.index_string(parameter.input_index ) +  ' - 1] + 1'
                    self.out + ')'
                else:    
                    self.out.n() + spec.input_var_name
                    self.out + '[' + self.index_string(parameter.input_index) + ']'
            if parameter.direction == RemoteFunction.INOUT:
                self.out.n() + '&' + spec.input_var_name 
                self.out + '[' + self.index_string(parameter.input_index) + ']'
            elif parameter.direction == RemoteFunction.OUT:
                self.out.n() + '&' + spec.output_var_name
                self.out + '[' + self.index_string(parameter.output_index) + ']'
    
        self.out.dedent()
        
    def output_lines_with_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
        
            if parameter.direction == RemoteFunction.INOUT:
                self.out.n() + spec.output_var_name
                self.out + '[' + self.index_string(parameter.output_index) + ']'
                self.out + ' = '
                self.out + spec.input_var_name + '[' + self.index_string(parameter.input_index) + ']'+';'
                
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
            self.out + 'reply_header.' + spec.counter_name 
            self.out + ' = ' + count + ';'
            pass
            
    def output_function_end(self):
        if len(self.specification.parameters) > 0:
            self.out.n()
            
        self.out + ')' + ';'
        
    def output_function_start(self):
        self.out.n() 
        if not self.specification.result_type is None:
            spec = self.dtype_to_spec[self.specification.result_type]
            self.out + spec.output_var_name
            self.out + '[' +  self.index_string(0) + ']' + ' = '
        self.out + self.specification.name + '('
        
    def output_casestmt_start(self):
        self.out + 'case ' + self.specification.id + ':'
        
    def output_casestmt_end(self):
        self.out.n() + 'break;'
        
        

class MakeACHeaderDefinitionStringOfALegacyFunctionSpecification(MakeCCodeString):

        
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
            self.out + ' '
            if parameter.is_output():
                self.out + '*' + ' '
            if parameter.datatype == 'string':
                self.out + '*' + ' '
            self.out + parameter.name
                
            
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
        self.out + self.specification.name + '('
        
class MakeACStringOfALegacyGlobalSpecification(MakeCCodeString):
    
    
            
    def start(self):
        self.output_casestmt_start()
        self.out.indent()
        
        spec = self.dtype_to_spec[self.legacy_global.datatype]
        self.out.n() + 'if(request_header.' + spec.counter_name
        self.out + ' == ' + 1 + '){'
        self.out.indent()
        self.out.n() + self.legacy_global.name + ' = ' 
        self.out + spec.input_var_name  + '[0]' + ';'
        self.out.dedent()
        self.out.n() + '} else {'
        self.out.indent()
        self.out.n() + 'reply_header.' + spec.counter_name
        self.out + ' = ' + 1 + ';'
        self.out.n() + spec.output_var_name + '[0]' 
        self.out + ' = ' + self.legacy_global.name + ';'
        self.out.dedent()
        self.out.n() + '}'
        
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
        
        
    def output_casestmt_start(self):
        self.out + 'case ' + self.legacy_global.id + ':'
        
    def output_casestmt_end(self):
        self.out.n() + 'break;'
        
        
        
class MakeCMessageHeaderClassDefinition(MakeCCodeString):
    @late
    def number_of_types(self):
        return len(self.dtype_to_spec)
        
    @late
    def length_of_the_header(self):
        return 2 + self.number_of_types
        
    def start(self):
        self.out + "class message_header {"
        self.out.lf() + "public:"
        self.out.indent()
        self.out.lf() + 'int tag;'
        self.out.lf() + 'int len;'
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'int ' + spec.counter_name;
            self.out + ';'
        
        self.make_constructor()
        self.make_send_function()
        self.make_recv_function()

        self.out.dedent()
        self.out.lf() + '}' + ';'
        
    def make_constructor(self):
        self.out.lf() + 'message_header(): tag(0), len(1)'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out + ', ' + spec.counter_name;
            self.out + '(0)'
        self.out + '{}'
        
    def make_destructor(self):
        self.out.lf() + '~message_header()'
        #for i, dtype in enumerate(dtypes):
        #    spec = self.dtype_to_spec[dtype]
        #    self.out + ', ' + spec.counter_name;
        #    self.out + '(0)'
        self.out + '{}'
        
    def make_send_function(self):
        self.out.lf() + 'void send(MPI::Intercomm & intercomm, int rank){'
        self.out.indent()
        self.out.lf() + 'int header[' + self.length_of_the_header + '];'
        self.out.lf() + 'header[0] = tag;'
        self.out.lf() + 'header[1] = len;'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype] 
            self.out.lf() + 'header[' + (i+2) + '] ='
            self.out + spec.counter_name + ';'
        self.out.lf() + 'intercomm.Send(header, '+ self.length_of_the_header 
        self.out +', MPI_INT, 0, 999);'
        self.out.dedent()
        self.out.lf() + '}'
        
    def make_recv_function(self):
        self.out.lf() + 'void recv(MPI::Intercomm & intercomm, int rank) {'
        self.out.indent()
        self.out.lf() + 'int header[' + self.length_of_the_header + '];'
        self.out.lf() + 'intercomm.Bcast(header, '+ self.length_of_the_header 
        self.out +', MPI_INT, 0);'
        self.out.lf() + 'tag = header[0];'
        self.out.lf() + 'len = header[1];'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype] 
            self.out.lf() + spec.counter_name + '='
            self.out + 'header[' + (i+2) + '];'
        self.out.dedent()
        self.out.lf() + '}'
        


class MakeACStringOfAClassWithLegacyFunctions\
    (MakeCodeStringOfAClassWithLegacyFunctions):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def make_legacy_function(self):
        return MakeACStringOfALegacyFunctionSpecification()
    
    def make_legacy_global(self):
        return MakeACStringOfALegacyGlobalSpecification()
        
    def start(self):
        self.output_mpi_include()
        self.output_local_includes()
        self.output_extra_content()
        
        self.out.lf().lf()
        self.output_header_class_definition()
        self.output_runloop_function_def_start()
        self.output_switch_start()
                
        self.output_legacy_functions()
        self.output_legacy_globals()
            
        self.output_switch_end()
        self.output_runloop_function_def_end()
        self.output_main()
        self._result = self.out.string
        
    def output_mpi_include(self):
        self.out.n() + '#include <mpi.h>'
        
    def output_local_includes(self):
        self.out.n()
        if hasattr(self.class_with_legacy_functions, 'include_headers'):
            for x in self.class_with_legacy_functions.include_headers:
                self.out.n() + '#include "' + x + '"'
    
    def output_extra_content(self):
        self.out.lf()
        if hasattr(self.class_with_legacy_functions, 'extra_content'):
            self.out.n() + self.class_with_legacy_functions.extra_content
    
            
    def output_header_class_definition(self):
        uc = MakeCMessageHeaderClassDefinition()
        uc.out = self.out
        uc.start()
        #self.out + HEADER_CLASS_STRING
        
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'void run_loop() {'
        self.out.indent()
        self.out.lf().lf() + 'MPI::Intercomm parent = '
        self.out + 'MPI::COMM_WORLD.Get_parent();'
        self.out.n() + 'int rank = parent.Get_rank();'
        self.out.lf().lf() + 'bool must_run_loop = true;'
        self.out.lf() + 'char * characters = 0;'
        
        self.out.lf().lf() + 'int max_len = ' + 10 + ';'
        self.output_new_statements(True)
        
        self.out.lf().lf() + 'while(must_run_loop) {'
        self.out.indent()
        
        
        self.out.lf()
        self.out.lf() + 'message_header request_header;'
        self.out.lf() + 'message_header reply_header;'
        self.out.lf()
        self.out.lf() + 'request_header.recv(parent,rank);'
            
        self.out.lf() + 'if (request_header.len > max_len) {'
        self.out.indent()
        self.out.lf() + 'max_len = request_header.len + 255;'
        self.output_delete_statements()
        self.output_new_statements(False)
        self.out.dedent()
        
        self.out.lf() + '}'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]    
            self.out.lf() 
            self.out + 'if(request_header.' + spec.counter_name + ' > 0) {'
            self.out.indent().lf() + 'parent.Bcast(' 
            self.out + spec.input_var_name 
            self.out + ', ' + 'request_header.' + spec.counter_name 
            self.out + ' * ' + 'request_header.len'
            self.out + ', ' + spec.mpi_type+ ', 0);'
            
            if dtype == 'string':
                self.out.lf() + 'characters = new char['
                self.out +  spec.input_var_name + '[' + 'request_header.' + spec.counter_name 
                self.out + ' * ' + 'request_header.len' + ' - 1] + 1'
                self.out + '];'
                self.out.lf() + 'parent.Bcast(' 
                self.out + 'characters'
                self.out + ',  ' 
                self.out + spec.input_var_name + '[' + 'request_header.' + spec.counter_name 
                self.out + ' * ' + 'request_header.len' + '- 1] + 1'
                self.out + ', ' + spec.mpi_type+ ', 0);'
            
            self.out.dedent().lf() +'}'
            
        self.out.lf().lf() + 'reply_header.tag = request_header.tag;'
        self.out.lf().lf() + 'reply_header.len = request_header.len;'
        
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
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]    
            self.out.lf() + 'if(reply_header.' 
            self.out + spec.counter_name + ' > 0) {'
            self.out.indent().lf() + 'parent.Send(' + spec.output_var_name 
            self.out + ', ' + 'reply_header.' + spec.counter_name 
            self.out + ' * ' + 'request_header.len'
            self.out + ', ' + spec.mpi_type+ ', 0, 999);'
            self.out.dedent().lf() +'}'
        
        self.out.lf() + 'if (characters) { delete characters; characters = 0;}'
        self.out.dedent()
        self.out.lf() + '}'
        
        self.output_delete_statements()
        
        self.out.lf().lf() + 'parent.Free();'
        self.out.dedent()
        self.out.lf() + '}'
    
    def output_delete_statements(self):
        for dtype_spec in self.dtype_to_spec.values():
            self.out.lf() + 'delete ' 
            self.out + dtype_spec.input_var_name  + ';'
            self.out.lf() + 'delete '
            self.out + dtype_spec.output_var_name  + ';'
            
    def output_new_statements(self, must_add_type):
        maximum_number_of_inputvariables_of_a_type = 255
        for dtype_spec in self.dtype_to_spec.values():
            self.out.lf() 
            if must_add_type:
                self.out + dtype_spec.type + ' * ' 
            self.out + dtype_spec.input_var_name 
            self.out + ' = new ' + dtype_spec.type 
            self.out + '[' + ' max_len * ' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
            self.out.lf() 
            if must_add_type:
                self.out + dtype_spec.type + ' * ' 
            self.out + dtype_spec.output_var_name 
            self.out + ' = new ' + dtype_spec.type 
            self.out + '[' + ' max_len * ' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
            
    def output_main(self):
        self.out.lf().lf() + 'int main(int argc, char *argv[])'
        self.out.lf() + '{'
        self.out.indent().lf() 
        self.out.lf() + 'int provided;'
        self.out.lf() + 'MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);'
        self.out.lf().lf() + 'run_loop();'
        self.out.lf().lf() + 'MPI_Finalize();'
        self.out.lf() + 'return 0;'
        self.out.dedent().lf()+'}'
    
        
        
        

class MakeACHeaderStringOfAClassWithLegacyFunctions\
    (MakeCodeStringOfAClassWithLegacyFunctions):

    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def make_legacy_function(self):
        return MakeACHeaderDefinitionStringOfALegacyFunctionSpecification()
        
    def start(self):  
        self.out + 'extern "C" {'
        self.out.indent().lf()
        self.output_legacy_functions()
        self._result = self.out
        self.out.dedent().lf() + '}'
    
        
        
        
