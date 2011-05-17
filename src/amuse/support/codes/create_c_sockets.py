from amuse.support.core import late
from amuse.support import exceptions
from amuse.support.codes.core import LegacyFunctionSpecification
from amuse.support.codes.create_code import GenerateASourcecodeString
from amuse.support.codes.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.support.codes.create_code import DTypeSpec, dtypes, DTypeToSpecDictionary

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('ints_in','ints_out',
                    'number_of_ints', 'int', 'MPI_INT'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'number_of_doubles', 'double', 'MPI_DOUBLE'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'number_of_floats', 'float', 'MPI_FLOAT'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'number_of_strings', 'int', 'MPI_INTEGER'),
    'bool' : DTypeSpec('bools_in', 'bools_out',
                    'number_of_bools', 'int', 'MPI_INTEGER'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'number_of_longs', 'long long int', 'MPI_LONG_LONG_INT'),
})

    
redirect_outputs_function_template = """
int internal__redirect_outputs(const char * stdoutfile, const char * stderrfile)
{
    char fullname[1024];
    int mpi_rank, mpi_err;
    
    fclose(stdin);
    
    mpi_err = MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    if (strcmp(stdoutfile, "none"))
    {
        if (strcmp(stdoutfile, "/dev/null") ) 
        {
            snprintf( fullname, 1024, "%s.%.3d",stdoutfile, mpi_rank );
        }
        else
        {
            snprintf( fullname, 1024, "%s",stdoutfile);
        }
        freopen(fullname, "a+", stdout);
    }
    
    if (strcmp(stderrfile, "none"))
    {         
        if (strcmp(stderrfile, "/dev/null") ) 
        {
            snprintf( fullname, 1024, "%s.%.3d",stderrfile, mpi_rank ); 
        }
        else
        {
            snprintf( fullname, 1024, "%s",stderrfile);
        }
        freopen(fullname, "a+", stderr);
    }
        
    return 0;
}
"""

EXIT_HANDLER_STRING = """\
void onexit(void) {
    int flag = 0;
    MPI_Finalized(&flag);
    
    if(!flag) {
        MPI::Intercomm parent = MPI::COMM_WORLD.Get_parent();
        int rank = parent.Get_rank();
        
        message_header reply_header;
        reply_header.tag = -2;
        reply_header.send(parent, rank);
        
        parent.Disconnect();
        MPI_Finalize();
    }
}
"""

class MakeCCodeString(GenerateASourcecodeString):
    @late
    def dtype_to_spec(self):
        return dtype_to_spec
       
         

class GenerateACStringOfAFunctionSpecification(MakeCCodeString):
    @late
    def specification(self):
        raise exceptions.AmuseException("No specification set, please set the specification first")
   
        
    def start(self):
        
        self.specification.prepare_output_parameters()
        self.output_casestmt_start()
        self.out.indent()
        
        if self.specification.must_handle_array:
            pass
        elif self.specification.can_handle_array:
            self.out.lf() + 'for (int i = 0 ; i < request_header.len; i++){'
            self.out.indent()
 
        
        self.output_lines_before_with_inout_variables()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'for (int i = 1 ; i < request_header.len; i++){'
                self.out.indent()
                self.out.lf() + spec.output_var_name + '[i]' + ' = ' + spec.output_var_name + '[0]' + ';'
                self.out.dedent()
                self.out.lf() + '}'
        elif self.specification.can_handle_array:
            self.out.dedent()
            self.out.lf() + '}'
        
        self.output_lines_with_number_of_outputs()
        self.output_casestmt_end()
        self.out.dedent()
        self._result = self.out.string
    
    def index_string(self, index, must_copy_in_to_out = False):
        if self.specification.must_handle_array and not must_copy_in_to_out:
            if index == 0:
                return '0'
            else:
                return '( %d * request_header.len)' % index
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                return '( %d * request_header.len) + i' % index
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
    
    def output_function_parameters(self):
        self.out.indent()
        
        first = True
        
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if first:
                first = False
            else:
                self.out + ' ,'
                
            if parameter.direction == LegacyFunctionSpecification.IN:
                if parameter.datatype == 'string':
                    self.out.n() + 'characters + ' 
                    self.out + '( ' + self.index_string(parameter.input_index) + '- 1 < 0 ? 0 :' 
                    self.out + spec.input_var_name
                    self.out + '[' + self.index_string(parameter.input_index ) +  ' - 1] + 1'
                    self.out + ')'
                else:    
                    self.input_var(spec.input_var_name, parameter.input_index)
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string':
                    self.output_var('output_strings', parameter.output_index)
                else:
                    self.input_var(spec.input_var_name, parameter.input_index)
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                if parameter.datatype == 'string': 
                    self.output_var('output_strings', parameter.output_index)
                else:
                    self.output_var(spec.output_var_name, parameter.output_index)
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out.n() + 'request_header.len'
    
        self.out.dedent()
    
    def output_lines_before_with_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
        
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string': 
                    self.out.n() + 'output_strings[' + self.index_string(parameter.output_index)  +'] = (characters + ' 
                    self.out + '( ' + self.index_string(parameter.input_index) + '- 1 < 0 ? 0 :' 
                    self.out + spec.input_var_name
                    self.out + '[' + self.index_string(parameter.input_index ) +  ' - 1] + 1'
                    self.out + '));'
    
    def output_lines_with_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'for (int i = 0 ; i < request_header.len; i++){'
                    self.out.indent()

                self.out.n() + spec.output_var_name
                self.out + '[' + self.index_string(parameter.output_index, must_copy_in_to_out = True) + ']'
                self.out + ' = '
                self.out + spec.input_var_name + '[' + self.index_string(parameter.input_index, must_copy_in_to_out = True) + ']'+';'
            
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
        
        
        
class MakeACStringOfALegacyGlobalSpecification(MakeCCodeString):
    @late
    def legacy_global(self):
        raise exceptions.AmuseException("No legacy_global set, please set the legacy_global first")
    
            
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
        
        
        
class GenerateCMessageHeaderClassDefinition(MakeCCodeString):
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
        


class GenerateACSourcecodeStringFromASpecificationClass\
    (GenerateASourcecodeStringFromASpecificationClass):

    @late
    def specification_class(self):
        raise exceptions.AmuseException("No specification_class set, please set the specification_class first")
    
    @late
    def dtype_to_spec(self):
        return dtype_to_spec

    def output_sourcecode_for_function(self):
        return GenerateACStringOfAFunctionSpecification()
    
    def make_legacy_global(self):
        return MakeACStringOfALegacyGlobalSpecification()
        
    def start(self):
        self.output_mpi_include()
        self.output_local_includes()
        self.output_extra_content()
        
        self.out.lf().lf()
        self.output_header_class_definition()
        self.output_onexit_fuinction()
        self.output_redirect_output()
        
        self.output_runloop_function_def_start()
        self.output_switch_start()
                
        self.output_sourcecode_for_functions()
        self.output_sourcecode_for_globals()
            
        self.output_switch_end()
        self.output_runloop_function_def_end()
        
        self.output_main()
        self._result = self.out.string
        
    def output_mpi_include(self):
        self.out.n() + '#include <mpi.h>'
        self.out.n() + '#include <iostream>'
        self.out.n() + '#include <string.h>'
        self.out.n() + '#include <stdlib.h>'
        self.out.n() + '#include <stdio.h>'
        
        
    def output_local_includes(self):
        self.out.n()
        if hasattr(self.specification_class, 'include_headers'):
            for x in self.specification_class.include_headers:
                self.out.n() + '#include "' + x + '"'
    
    def output_extra_content(self):
        self.out.lf()
        if hasattr(self.specification_class, 'extra_content'):
            self.out.n() + self.specification_class.extra_content
    
            
    def output_header_class_definition(self):
        uc = GenerateCMessageHeaderClassDefinition()
        uc.out = self.out
        uc.start()
    
    def output_onexit_fuinction(self):
        self.out.lf() + EXIT_HANDLER_STRING
        
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'void run_loop() {'
        self.out.indent()
        self.out.lf().lf() + 'MPI::Intercomm parent = '
        self.out + 'MPI::COMM_WORLD.Get_parent();'
        self.out.n() + 'int rank = parent.Get_rank();'
        self.out.lf().lf() + 'bool must_run_loop = true;'
        self.out.lf() + 'char * characters = 0, * output_characters = 0;'
        
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
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)
            if max == 0:
                continue
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
                self.out + ', ' +'MPI_CHARACTER'+ ', 0);'
            
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
        self.out.lf().lf() + 'MPI::COMM_WORLD.Barrier();'
        self.out.lf().lf()
        
        self.out.lf() + 'if(rank == 0) {'
        
        self.out.indent().lf()
        self.out.lf().lf() + 'reply_header.send(parent, rank);'
        self.out.lf().lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]  
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            if max == 0:
                continue
                  
            self.out.lf() + 'if(reply_header.' 
            self.out + spec.counter_name + ' > 0) {'
            self.out.indent().lf()
            if dtype == 'string':
                self.out.lf() + 'int offset = 0;'
                self.out.lf() + 'for( int i = 0; i < ' +  '(reply_header.' + spec.counter_name + '*' + 'request_header.len) ; i++) {'
                self.out.indent().lf()
                self.out.lf() + 'int length = strlen(' + 'output_strings' + '[i]);'
                self.out.lf() + 'offset += length;'
                self.out.lf() + spec.output_var_name + '[i] = offset;'
                self.out.lf() + 'offset += 1;'
                self.out.dedent().lf() +'}'
                self.out.lf() + 'output_characters = new char[offset + 1];'
                
                self.out.lf() + 'offset = 0;'
                self.out.lf() + 'for( int i = 0; i < ' + '(reply_header.' + spec.counter_name + '*' + 'request_header.len)  ; i++) {'
                self.out.indent().lf()
                self.out.lf() + 'strcpy(output_characters+offset, output_strings[i]);'
                self.out.lf() + 'offset = ' + spec.output_var_name + '[i] + 1;'
                self.out.dedent().lf() +'}'
                self.out.lf()
                
            self.out + 'parent.Send(' + spec.output_var_name 
            self.out + ', ' + 'reply_header.' + spec.counter_name 
            self.out + ' * ' + 'request_header.len'
            self.out + ', ' + spec.mpi_type+ ', 0, 999);'
            
            if dtype == 'string':
                self.out.lf() + 'parent.Send(' + 'output_characters'
                self.out + ',  offset' # + spec.output_var_name +'['+'reply_header.' + spec.counter_name +' - 1] + 1'
                self.out + ', ' + 'MPI_BYTE' + ', 0, 999);'
                
            self.out.dedent().lf() +'}'
            
        self.out.dedent().lf()
        self.out.lf() + '}'
        
        self.out.lf() + 'if (characters) { delete[] characters; characters = 0;}'
        self.out.lf() + 'if (output_characters) { delete[] output_characters; output_characters = 0;}'
        self.out.dedent()
        self.out.lf() + '}'
        
        self.output_delete_statements()
        
        self.out.lf().lf() + 'parent.Disconnect();'
        self.out.dedent()
        self.out.lf() + '}'
    
    def output_delete_statements(self):
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
            if max > 0:
                self.out.lf() + 'delete[] ' 
                self.out + dtype_spec.input_var_name  + ';'
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            if max > 0:
                self.out.lf() + 'delete[] '
                self.out + dtype_spec.output_var_name  + ';'
        
        max = 0
        dtype='string'
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
        if max > 0:
            self.out.lf() + 'delete[] '
            self.out + 'output_strings' + ';'
            
    def output_new_statements(self, must_add_type):
        maximum_number_of_inputvariables_of_a_type = 255
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            self.out.lf() 
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype,0)
            if max > 0:
                self.output_new_statement(
                    must_add_type,
                    dtype_spec.type,
                    dtype_spec.input_var_name,
                    max
                )
                
            max = 0
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
            if max > 0:
                self.output_new_statement(
                    must_add_type,
                    dtype_spec.type,
                    dtype_spec.output_var_name,
                    max
                )
                
        max = 0
        dtype = 'string'
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype,0)
        if max > 0:
            self.output_new_statement(
                must_add_type,
                'char *',
                'output_strings',
                max
            )
            
    def output_new_statement(self, must_add_type,  type, var_name, maximum_number_of_inputvariables_of_a_type):
        if must_add_type:
            self.out + type + ' * ' 
            
        self.out + var_name
        self.out + ' = new ' + type
        self.out + '[' + ' max_len * ' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
        self.out.lf() 
        
    def output_redirect_output(self):
        self.out.lf().lf() + redirect_outputs_function_template
    
    def output_main(self):
        self.out.lf().lf() + 'int main(int argc, char *argv[])'
        self.out.lf() + '{'
        self.out.indent().lf() 
        #self.out.lf() + 'int provided;'
        self.out.lf() + 'MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);'
        #self.out.lf() + 'MPI::Init(argc, argv);'
        self.out.lf() + 'atexit(onexit);'
        self.out.lf().lf() + 'run_loop();'
        self.out.lf().lf() + 'MPI_Finalize();'
        self.out.lf() + 'return 0;'
        self.out.dedent().lf()+'}'
    
        
        
        

        
        
        
