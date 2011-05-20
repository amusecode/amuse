from amuse.support.core import late
from amuse.support import exceptions
from amuse.support.codes.core import LegacyFunctionSpecification
from amuse.support.codes.create_code import GenerateASourcecodeString
from amuse.support.codes.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.support.codes.create_code import DTypeSpec, DTypeToSpecDictionary

dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('ints_in', 'ints_out',
                    'number_of_ints', 'int', 'MPI_INT'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'number_of_longs', 'long long int', 'MPI_LONG_LONG_INT'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'number_of_floats', 'float', 'MPI_FLOAT'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'number_of_doubles', 'double', 'MPI_DOUBLE'),
    'bool' : DTypeSpec('booleans_in', 'booleans_out',
                    'number_of_booleans', 'int', 'MPI_INTEGER'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'number_of_strings', 'int', 'MPI_INTEGER'),
})

dtypes = ['int32', 'int64', 'float32', 'float64', 'bool', 'string']

REDIRECT_OUTPUTS_FUNCTION_STRING = """
int internal__redirect_outputs(const char * stdoutfile, const char * stderrfile)
{
        
    return 0;
}
"""

EXIT_HANDLER_FUNCTION_STRING = """\
void onexit(void) {

    close(socketfd);
        
}
"""

SEND_FUNCTION_STRING = """\
void send(void *buffer, int length, int file_descriptor) {
    int total_written = 0;
    int written;

    while (total_written < length) {
        written = write(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written);

        if (written == -1) {
            fprintf(stderr, "could not write data\\n");
            exit(1);
        }

        total_written = total_written + written;
    }
}
"""

RECEIVE_FUNCTION_STRING = """\
void receive(void *buffer, int length, int file_descriptor) {
    int total_read = 0;
    int bytes_read;

    while (total_read < length) {
        bytes_read = read(file_descriptor, ((char *) buffer) + total_read,
                        length - total_read);

        if (bytes_read == -1) {
            fprintf(stderr, "could not read data\\n");
            exit(1);
        }

        total_read = total_read + bytes_read;
    }
}
"""

MAIN_FUNCTION_STRING = """\
int main(int argc, char *argv[]) {
    int portno;
    struct sockaddr_in serv_addr;
    struct hostent *server;

    fprintf(stderr, "main!\\n");

    portno = atoi(argv[1]);
    socketfd = socket(AF_INET, SOCK_STREAM, 0);

    if (socketfd < 0) {
            fprintf(stderr, "cannot open socket\\n");
            exit(0);
    }

    server = gethostbyname("localhost");

    fprintf(stderr, "connecting...\\n");

    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr,
                    server->h_length);
    serv_addr.sin_port = htons(portno);
    if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
                    < 0) {
            fprintf(stderr, "cannot connect socket\\n");
            exit(0);

    }

    fprintf(stderr, "running...\\n");

    run_loop(socketfd);

    close(socketfd);

    return 0;
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
            self.out.lf() + 'for (int i = 0 ; i < request_header.length; i++){'
            self.out.indent()
 
        
        self.output_lines_before_with_inout_variables()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_lines_with_inout_variables()
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'for (int i = 1 ; i < request_header.length; i++){'
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
    
    def index_string(self, index, must_copy_in_to_out=False):
        if self.specification.must_handle_array and not must_copy_in_to_out:
            if index == 0:
                return '0'
            else:
                return '( %d * request_header.length)' % index
        elif self.specification.can_handle_array or (self.specification.must_handle_array and must_copy_in_to_out):
            if index == 0:
                return 'i'
            else:
                return '( %d * request_header.length) + i' % index
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
                    self.input_var('characters_in', parameter.input_index)
                else:    
                    self.input_var(spec.input_var_name, parameter.input_index)
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string':
                    self.output_var('characters_out', parameter.output_index)
                else:
                    self.input_var(spec.input_var_name, parameter.input_index)
            elif parameter.direction == LegacyFunctionSpecification.OUT:
                if parameter.datatype == 'string': 
                    self.output_var('characters_out', parameter.output_index)
                else:
                    self.output_var(spec.output_var_name, parameter.output_index)
            elif parameter.direction == LegacyFunctionSpecification.LENGTH:
                self.out.n() + 'request_header.length'
    
        self.out.dedent()
    
    def output_lines_before_with_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
        
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string': 
                    self.out.n() + 'output_strings[' + self.index_string(parameter.output_index) + '] = (characters + ' 
                    self.out + '( ' + self.index_string(parameter.input_index) + '- 1 < 0 ? 0 :' 
                    self.out + spec.input_var_name
                    self.out + '[' + self.index_string(parameter.input_index) + ' - 1] + 1'
                    self.out + '));'
    
    def output_lines_with_inout_variables(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'for (int i = 0 ; i < request_header.length; i++){'
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
            self.out + 'reply_header.' + spec.counter_name 
            self.out + ' = ' + count + ' * request_header.length;'
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
            self.out + '[' + self.index_string(0) + ']' + ' = '
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
        self.out + spec.input_var_name + '[0]' + ';'
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
        return 4 + self.number_of_types
        
    def start(self):
        self.out + 'static int HEADER_SIZE = ' + self.length_of_the_header + ' * sizeof(int);'
        self.out.lf() 
        
        self.out.lf() + "class message_header {"
        self.out.lf()
        self.out.indent()
        self.out.lf() + "public:"
        self.out.lf() + 'int flags;'
        self.out.lf() + 'int id;'
        self.out.lf() + 'int function_id;'
        self.out.lf() + 'int length;'
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'int ' + spec.counter_name;
            self.out + ';'

        self.out.lf()
        
        self.make_constructor()
        self.out.lf()

        self.make_send_function()
        self.out.lf()

        self.make_recv_function()
        self.out.lf()

        self.out.dedent()
        self.out.lf() + '};'
        self.out.lf()
        
    def make_constructor(self):
        self.out.lf() + 'message_header(): flags(0), id(0), function_id(0), length(1)'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out + ', ' + spec.counter_name;
            self.out + '(0)'
        self.out + ' {'
        self.out.lf() + '}'

        
    def make_destructor(self):
        self.out.lf() + '~message_header()'
        #for i, dtype in enumerate(dtypes):
        #    spec = self.dtype_to_spec[dtype]
        #    self.out + ', ' + spec.counter_name;
        #    self.out + '(0)'
        self.out + '{}'
        
    def make_send_function(self):
        self.out.lf() + 'void send(int file_descriptor) {'
        self.out.indent()
        self.out.lf() + 'int header[' + self.length_of_the_header + '];'
        self.out.lf() + 'int total_written;'
        self.out.lf() + 'int bytes_written;'
        self.out.lf()
        
        self.out.lf() + 'header[0] = flags;'
        self.out.lf() + 'header[1] = id;'
        self.out.lf() + 'header[2] = function_id;'
        self.out.lf() + 'header[3] = length;'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype] 
            self.out.lf() + 'header[' + (i + 4) + '] = '
            self.out + spec.counter_name + ';'


        self.out.lf()
        self.out.lf() + 'total_written = 0;'

        self.out.lf() + 'while (total_written < HEADER_SIZE) {'
        self.out.indent();

        self.out.lf() + 'bytes_written = write(file_descriptor, &header + total_written,'
        self.out.lf() + '                      HEADER_SIZE - total_written);'
        
        self.out.lf() + 'if (bytes_written == -1) {'
        self.out.indent()
        self.out.lf() + 'fprintf(stderr, "could not write header\\n");'
        self.out.lf() + 'exit(1);'
        self.out.dedent()
        self.out.lf() + '}'

        self.out.lf() + 'total_written = total_written + bytes_written;'
                        
        self.out.dedent()
        self.out.lf() + '}'
        self.out.dedent()
        self.out.lf() + '}'            
            
    def make_recv_function(self):
        self.out.lf() + 'void recv(int file_descriptor) {'
        self.out.indent()
        self.out.lf() + 'int header[' + self.length_of_the_header + '];'
        self.out.lf() + 'int total_read;'
        self.out.lf() + 'int bytes_read;'
        
        self.out.lf()
        self.out.lf() + 'total_read = 0;'

        self.out.lf() + 'while (total_read < HEADER_SIZE) {'
        self.out.indent();

        self.out.lf() + 'bytes_read = read(file_descriptor, &header + total_read,'
        self.out.lf() + '                  HEADER_SIZE - total_read);'
        
        self.out.lf() + 'if (bytes_read == -1) {'
        self.out.indent()
        self.out.lf() + 'fprintf(stderr, "could not read header\\n");'
        self.out.lf() + 'exit(1);'
        self.out.dedent()
        self.out.lf() + '}'

        self.out.lf() + 'total_read = total_read + bytes_read;'
                        
        self.out.dedent()
        self.out.lf() + '}'            

        self.out.lf()
        self.out.lf() + 'flags = header[0];'
        self.out.lf() + 'id = header[1];'
        self.out.lf() + 'function_id = header[2];'
        self.out.lf() + 'length = header[3];'
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype] 
            self.out.lf() + spec.counter_name + ' = '
            self.out + 'header[' + (i + 4) + '];'
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
        self.output_includes()
        self.output_local_includes()
        self.output_global_variables()
        self.output_extra_content()
        
        self.output_header_class_definition()
        self.output_redirect_outputs_function()
        self.output_onexit_function()
        self.output_send_function()
        self.output_receive_function()
        
        #start of main loop
        
        self.output_runloop_function_def_start()
        
        self.output_switch_start()
                
        self.output_sourcecode_for_functions()
        self.output_sourcecode_for_globals()
            
        self.output_switch_end()
        self.output_runloop_function_def_end()
        
        #end of main loop
        
        self.output_main()
        self._result = self.out.string
        
    def output_global_variables(self):
        self.out.lf()
        self.out.lf() + 'static int ERROR_FLAG = 256;'
        self.out.lf() + 'int socketfd = 0;' 
        
    def output_includes(self):
        self.out.n() + '#include <mpi.h>'
        self.out.n() + '#include <iostream>'
        self.out.n() + '#include <string.h>'
        self.out.n() + '#include <stdlib.h>'
        self.out.n() + '#include <stdio.h>'
        self.out.n() + '#include <sys/socket.h>'
        self.out.n() + '#include <netinet/in.h>'
        self.out.n() + '#include <netdb.h>'
        
        
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
        
    def output_redirect_outputs_function(self):
        self.out.lf() + REDIRECT_OUTPUTS_FUNCTION_STRING
    
    def output_onexit_function(self):
        self.out.lf() + EXIT_HANDLER_FUNCTION_STRING
        
    def output_send_function(self):
        self.out.lf() + SEND_FUNCTION_STRING
        
    def output_receive_function(self):
        self.out.lf() + RECEIVE_FUNCTION_STRING
        
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'void run_loop(int socketfd) {'
        self.out.indent()
        self.out.lf().lf() + 'bool must_run_loop = true;'
        
        self.out.lf().lf() + 'int max_len = 10;'
        self.out.lf();
        self.output_new_statements(True)
        
        self.out.lf().lf() + 'while(must_run_loop) {'
        self.out.indent()
        
        self.out.lf()
        self.out.lf() + 'message_header request_header;'
        self.out.lf() + 'message_header reply_header;'
        self.out.lf()
        self.out.lf() + 'request_header.recv(socketfd);'
        self.out.lf()
            
        self.out.lf() + 'if (request_header.length > max_len) {'
        self.out.indent()
        self.out.lf() + 'max_len = request_header.length + 255;'
        self.output_delete_statements()
        self.output_new_statements(False)
        self.out.dedent()
        
        self.out.lf() + '}'
        self.out.lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]    
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
            if max == 0:
                continue
            self.out.lf() 
            self.out + 'if(request_header.' + spec.counter_name + ' > 0) {'
            self.out.indent().lf() + 'receive(' 
            self.out + spec.input_var_name 
            self.out + ', ' + 'request_header.' + spec.counter_name 
            self.out + ' * sizeof(' + spec.type + '), socketfd);'
            
            if dtype == 'string':
                self.out.lf() + 'for (int i = 0; i < request_header.number_of_strings; i++) {'
                self.out.indent()
                self.out.lf() + 'characters_in[i] = new char[strings_in[i]];'
                self.out.lf() + 'receive(characters_in[i], strings_in[i], socketfd);'
                self.out.dedent().lf() + '}'
                                      
            self.out.dedent().lf() + '}'
            self.out.lf() 
            
        self.out.lf() + 'reply_header.function_id = request_header.function_id;'
        self.out.lf() + 'reply_header.length = request_header.length;'
        
    def output_switch_start(self):
        self.out.lf().lf() + 'switch(request_header.function_id) {'
        self.out.indent()
        self.out.lf() + 'case 0:'
        self.out.indent().lf() + 'must_run_loop = false;'
        self.out.lf() + 'break;'
        self.out.dedent()
        
    def output_switch_end(self):
        self.out.lf() + 'default:'
        self.out.indent().lf() + 'reply_header.flags = reply_header.flags & ERROR_FLAG;'
        self.out.dedent()
        self.out.dedent().lf() + '}'
        
    def output_runloop_function_def_end(self):
        self.out.lf()
        self.out.lf().lf() + 'reply_header.send(socketfd);'
        self.out.lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]  
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            if max == 0:
                continue
                  
            self.out.lf() + 'if(reply_header.' 
            self.out + spec.counter_name + ' > 0) {'
            self.out.indent().lf()
                
            self.out + 'send(' + spec.output_var_name + ', ' + 'reply_header.' + spec.counter_name + ' * sizeof(' + spec.type + '), socketfd);'
            
            self.out.dedent().lf() + '}'
            self.out.lf()
            
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get('string', 0)
        if max != 0:
            self.out.lf() + 'for (int i = 0; i < reply_header.number_of_strings; i++) {'
            self.out.indent()
            self.out.lf() + 'send(characters_out[i], reply_header.strings_out[i] * sizeof(char), socketfd);'
            self.out.dedent();
            self.out.lf() + '}'
            self.out.lf()

            
        max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get('string', 0)
        if max != 0:
            self.out.lf() + 'for (int i = 0; i < request_header.number_of_strings; i++) {'
            self.out.indent()
            self.out.lf() + 'delete[] characters_in[i];'
            self.out.dedent();
            self.out.lf() + '}'
            self.out.lf()
            
        self.out.dedent()
        self.out.lf() + '}'
        self.output_delete_statements()
        
        self.out.dedent()
        self.out.lf() + '}'
        self.out.lf()
    
    def output_delete_statements(self):
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
            if max > 0:
                self.out.lf() + 'delete[] ' 
                self.out + dtype_spec.input_var_name + ';'
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            if max > 0:
                self.out.lf() + 'delete[] '
                self.out + dtype_spec.output_var_name + ';'
        
        max = 0
        dtype = 'string'
        max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
        if max > 0:
            self.out.lf() + 'delete[] '
            self.out + 'characters_in' + ';'
            
        max = 0
        dtype = 'string'
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
        if max > 0:
            self.out.lf() + 'delete[] '
            self.out + 'characters_out' + ';'
            
            
    def output_new_statements(self, must_add_type):
        maximum_number_of_inputvariables_of_a_type = 255
        
        for dtype in self.dtype_to_spec.keys():
            dtype_spec = self.dtype_to_spec[dtype]
            self.out.lf() 
            max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
            if max > 0:
                self.output_new_statement(
                    must_add_type,
                    dtype_spec.type,
                    dtype_spec.input_var_name,
                    max
                )
                
            max = 0
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            if max > 0:
                self.output_new_statement(
                    must_add_type,
                    dtype_spec.type,
                    dtype_spec.output_var_name,
                    max
                )

        #add space for characters too
        
        dtype = 'string'
        max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get(dtype, 0)
        if max > 0:
            self.output_new_statement(
                must_add_type,
                'char *',
                'characters_in',
                max
            )
                
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
        if max > 0:
            self.output_new_statement(
                must_add_type,
                'char *',
                'characters_out',
                max
            )
            
    def output_new_statement(self, must_add_type, type, var_name, maximum_number_of_inputvariables_of_a_type):
        if must_add_type:
            self.out + type + ' * ' 
            
        self.out + var_name
        self.out + ' = new ' + type
        self.out + '[' + ' max_len * ' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
        self.out.lf() 
    
    def output_main(self):
        self.out.lf() + MAIN_FUNCTION_STRING
    
        
        
        

        
        
        
