from amuse.support.core import late
from amuse.support import exceptions
from amuse.rfi.tools.create_code import GenerateASourcecodeString
from amuse.rfi.tools.create_code import GenerateASourcecodeStringFromASpecificationClass
from amuse.rfi.tools.create_code import DTypeSpec
from amuse.rfi.tools.create_code import DTypeToSpecDictionary
from amuse.rfi.core import LegacyFunctionSpecification
dtype_to_spec = DTypeToSpecDictionary({
    'int32' : DTypeSpec('ints_in', 'ints_out',
                    'HEADER_INTEGER_COUNT', 'int', 'MPI_INTEGER'),
    'int64' : DTypeSpec('longs_in', 'longs_out',
                    'HEADER_LONG_COUNT', 'long long int', 'MPI_INTEGER'),
    'float32' : DTypeSpec('floats_in', 'floats_out',
                    'HEADER_FLOAT_COUNT', 'float', 'MPI_FLOAT'),
    'float64' : DTypeSpec('doubles_in', 'doubles_out',
                    'HEADER_DOUBLE_COUNT', 'double', 'MPI_DOUBLE'),
    'bool' : DTypeSpec('booleans_in', 'booleans_out',
                    'HEADER_BOOLEAN_COUNT', 'int', 'MPI_INTEGER'),
    'string' : DTypeSpec('strings_in', 'strings_out',
                    'HEADER_STRING_COUNT', 'int', 'MPI_CHARACTER'),
})

dtypes = ['int32', 'int64', 'float32', 'float64', 'bool', 'string']

GLOBAL_VARIABLES_STRING = """
static int ERROR_FLAG = 256;
static int HEADER_SIZE = 10; //integers

static int HEADER_FLAGS = 0;
static int HEADER_CALL_ID = 1;
static int HEADER_FUNCTION_ID = 2;
static int HEADER_CALL_COUNT = 3;
static int HEADER_INTEGER_COUNT = 4;
static int HEADER_LONG_COUNT = 5;
static int HEADER_FLOAT_COUNT = 6;
static int HEADER_DOUBLE_COUNT = 7;
static int HEADER_BOOLEAN_COUNT = 8;
static int HEADER_STRING_COUNT = 9;

static bool TRUE_BYTE = 1;
static bool FALSE_BYTE = 0;

int socketfd = 0;
"""    

REDIRECT_OUTPUTS_FUNCTION_STRING = """
int internal__redirect_outputs(const char * stdoutfile, const char * stderrfile) {
    return 0;
}
"""

EXIT_HANDLER_FUNCTION_STRING = """\
void onexit(void) {
    close(socketfd);
}
"""

SEND_FUNCTION_STRING = """\
void send(void *buffer, int length, int file_descriptor, int rank) {
    int total_written = 0;
    int bytes_written;

    if (rank != 0) {
        return;
    }

    while (total_written < length) {
        bytes_written = write(file_descriptor, ((char *) buffer) + total_written,
                        length - total_written);

        if (bytes_written == -1) {
            fprintf(stderr, "could not write data\\n");
            exit(1);
        }

        total_written = total_written + bytes_written;
    }
}
"""

RECEIVE_FUNCTION_STRING = """\
void receive(void *buffer, int length, int file_descriptor, int rank) {
    int total_read = 0;
    int bytes_read;

    if (rank != 0) {
        return;
    }

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
    int rank;
    struct sockaddr_in serv_addr;
    struct hostent *server;
    
    MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
    
    rank = MPI::COMM_WORLD.Get_rank();

    if (rank == 0) {

        portno = atoi(argv[1]);
        socketfd = socket(AF_INET, SOCK_STREAM, 0);

        if (socketfd < 0) {
                fprintf(stderr, "cannot open socket\\n");
                exit(1);
        }

        server = gethostbyname("localhost");

        bzero((char *) &serv_addr, sizeof(serv_addr));
        serv_addr.sin_family = AF_INET;
        bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr,
                                    server->h_length);
        serv_addr.sin_port = htons(portno);
        if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
                    < 0) {
                fprintf(stderr, "cannot connect socket\\n");
                exit(1);

        }

        fprintf(stderr, "sockets: finished initializing code\\n");
        
    }

    run_loop(socketfd, rank);
    
    if (rank == 0) {
        close(socketfd);
    }

    MPI_Finalize();

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
            self.out.lf() + 'for (int i = 0 ; i < call_count; i++) {'
            self.out.indent()
 
        self.output_copy_inout_variables_before_function()
        self.output_function_start()
        self.output_function_parameters()
        self.output_function_end()
        self.output_copy_inout_variables_after_function()
        
        if self.specification.must_handle_array:
            if not self.specification.result_type is None:
                spec = self.dtype_to_spec[self.specification.result_type]
                self.out.lf() + 'for (int i = 1 ; i < call_count; i++) {'
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
                self.out.n() + 'call_count'
    
        self.out.dedent()
    
    
        
    def output_copy_inout_variables_before_function(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
        
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if parameter.datatype == 'string': 
                    self.out.n() + 'characters_out[' + self.index_string(parameter.output_index, must_copy_in_to_out=True) + ']'
                    self.out + ' = '
                    self.out + 'characters_in[' + self.index_string(parameter.input_index, must_copy_in_to_out=True) + ']'
                    self.out + ';'

    def output_copy_inout_variables_after_function(self):
        for parameter in self.specification.parameters:
            spec = self.dtype_to_spec[parameter.datatype]
            
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                if self.specification.must_handle_array:
                    self.out.lf() + 'for (int i = 0 ; i < call_count; i++) {'
                    self.out.indent()

                if parameter.datatype != 'string':
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
        self.out.n() + 'if(header_in[' + spec.counter_name
        self.out + '] == ' + 1 + '){'
        self.out.indent()
        self.out.n() + self.legacy_global.name + ' = ' 
        self.out + spec.input_var_name + '[0]' + ';'
        self.out.dedent()
        self.out.n() + '} else {'
        self.out.indent()
        self.out.n() + 'header_out[' + spec.counter_name
        self.out + '] = ' + 1 + ';'
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
        self.out.lf() + GLOBAL_VARIABLES_STRING;
        
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
        
    def output_redirect_outputs_function(self):
        self.out.lf() + REDIRECT_OUTPUTS_FUNCTION_STRING
    
    def output_onexit_function(self):
        self.out.lf() + EXIT_HANDLER_FUNCTION_STRING
        
    def output_send_function(self):
        self.out.lf() + SEND_FUNCTION_STRING
        
    def output_receive_function(self):
        self.out.lf() + RECEIVE_FUNCTION_STRING
        
    def output_runloop_function_def_start(self):
        self.out.lf().lf() + 'void run_loop(int socketfd, int rank) {'
        self.out.indent()
        self.out.lf().lf() + 'bool must_run_loop = true;'
        
        self.out.lf().lf() + 'int max_call_count = 10;'
        self.out.lf().lf() + 'int call_count;'
        self.out.lf() + 'char * buffer = new char[50];'

        self.out.lf();
        self.output_new_statements(True)
        
        self.out.lf()
        self.out.lf() + 'int header_in[HEADER_SIZE];'
        self.out.lf() + 'int header_out[HEADER_SIZE];'
        
        self.out.lf().lf() + 'while(must_run_loop) {'
        self.out.indent()

        self.out.lf() + 'receive(header_in, HEADER_SIZE * sizeof(int), socketfd, rank);'
        self.out.lf() + 'MPI::COMM_WORLD.Bcast(header_in, HEADER_SIZE, MPI_INT, 0);'
        
        self.out.lf() + '/*fprintf(stderr, "got header %d %d %d %d %d %d %d %d %d %d\\n", header_in[0], header_in[1], header_in[2], header_in[3], header_in[4], header_in[5], header_in[6], header_in[7], header_in[8], header_in[9]);*/'
        
        self.out.lf() + 'call_count = header_in[HEADER_CALL_COUNT];'
            
        self.out.lf() + 'if (call_count > max_call_count) {'
        self.out.indent()
        self.output_delete_statements()
        self.out.lf() + 'max_call_count = call_count + 255;'
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
            self.out + 'if(header_in[' + spec.counter_name + '] > 0) {'
            self.out.indent()
            
            
            if dtype == 'bool':
                self.out.lf() + 'for (int i = 0; i < header_in[HEADER_BOOLEAN_COUNT]; i++) {'
                self.out.indent()
                self.out.lf() + 'booleans_in[i] = 0;'
                self.out.lf() + 'receive(&booleans_in[i], 1, socketfd , rank);'
                self.out.dedent().lf() + '}'
                self.out.lf() + 'MPI::COMM_WORLD.Bcast(' + spec.input_var_name + ', header_in[' + spec.counter_name + '], ' + spec.mpi_type + ', 0);'

            elif dtype == 'string':
                self.out.lf() + 'receive(' 
                self.out + spec.input_var_name 
                self.out + ', ' + 'header_in[' + spec.counter_name 
                self.out + '] * sizeof(' + spec.type + '), socketfd, rank);'
                self.out.lf() + 'MPI::COMM_WORLD.Bcast(strings_in, header_in[HEADER_STRING_COUNT], MPI_INT, 0);'
                self.out.lf() + 'for (int i = 0; i < header_in[HEADER_STRING_COUNT]; i++) {'
                self.out.indent()
                self.out.lf() + 'characters_in[i] = new char[strings_in[i] + 1];'
                self.out.lf() + 'receive(characters_in[i], strings_in[i], socketfd, rank);'
                self.out.lf() + "characters_in[i][strings_in[i]] = '\\0';"
                self.out.lf() + "MPI::COMM_WORLD.Bcast(characters_in[i], strings_in[i], MPI_CHARACTER, 0);"
                self.out.dedent().lf() + '}'
            else:
                self.out.lf() + 'receive(' 
                self.out + spec.input_var_name 
                self.out + ', ' + 'header_in[' + spec.counter_name 
                self.out + '] * sizeof(' + spec.type + '), socketfd, rank);'
                self.out.lf() + 'MPI::COMM_WORLD.Bcast(' + spec.input_var_name + ', header_in[' + spec.counter_name + '], ' + spec.mpi_type + ', 0);'
            self.out.dedent().lf() + '}'
            self.out.lf() 

        self.out.lf() + 'header_out[HEADER_FLAGS] = 0;'
        self.out.lf() + 'header_out[HEADER_CALL_ID] = header_in[HEADER_CALL_ID];'
        self.out.lf() + 'header_out[HEADER_FUNCTION_ID] = header_in[HEADER_FUNCTION_ID];'
        self.out.lf() + 'header_out[HEADER_CALL_COUNT] = header_in[HEADER_CALL_COUNT];'
            
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]
            self.out.lf() + 'header_out[' + spec.counter_name + '] = 0;'
            
            
        self.out.lf() + 'MPI::COMM_WORLD.Barrier();'
        
    def output_switch_start(self):
        self.out.lf().lf() + 'switch(header_in[HEADER_FUNCTION_ID]) {'
        self.out.indent()
        self.out.lf() + 'case 0:'
        self.out.indent().lf() + 'must_run_loop = false;'
        self.out.lf() + 'break;'
        self.out.dedent()
        
    def output_switch_end(self):
        self.out.lf() + 'default:'
        self.out.indent().lf() + 'header_out[HEADER_FLAGS] = header_out[HEADER_FLAGS] | ERROR_FLAG;'
        self.out.lf() + 'characters_out[0] = buffer;'
        self.out.lf() + 'sprintf(buffer, "unknown function id: %d\\n", header_in[HEADER_FUNCTION_ID]);'
        self.out.lf() + 'fprintf(stderr, "unknown function id: %d\\n", header_in[HEADER_FUNCTION_ID]);'

        self.out.lf() + 'header_out[HEADER_STRING_COUNT] = 1;'
        
        self.out.dedent()
        self.out.dedent().lf() + '}'
        
    def output_runloop_function_def_end(self):
        self.out.lf() + 'MPI::COMM_WORLD.Barrier();'
        
        self.out.lf() + 'if (rank == 0) {'
        self.out.lf().indent()

        self.out.lf().lf() + 'send(header_out, HEADER_SIZE * sizeof(int), socketfd, rank);'
        self.out.lf()
        
        self.out.lf() + 'for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {'
        self.out.indent()
        self.out.lf() + 'strings_out[i] = strlen(characters_out[i]);'
        self.out.dedent();
        self.out.lf() + '}'
        self.out.lf()
        
        for i, dtype in enumerate(dtypes):
            spec = self.dtype_to_spec[dtype]  
            
            max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
            if max == 0:
                continue
                  
            self.out.lf() + 'if(header_out[' 
            self.out + spec.counter_name + '] > 0) {'
            self.out.indent().lf()
                
            if dtype == 'string':
                self.out + 'send(' + spec.output_var_name + ', ' + 'header_out[' + spec.counter_name + '] * sizeof(' + spec.type + '), socketfd, rank);'
                self.out.lf() + 'for (int i = 0; i < header_out[HEADER_STRING_COUNT]; i++) {'
                self.out.indent()
                self.out.lf() + 'send(characters_out[i], strings_out[i] * sizeof(char), socketfd, rank);'
                self.out.dedent();
                self.out.lf() + '}'
                self.out.lf()
                
            elif dtype == 'bool':
                self.out.lf() + 'for (int i = 0; i < header_out[HEADER_BOOLEAN_COUNT]; i++) {'
                self.out.indent()
                self.out.lf() + 'if (booleans_out[i]) {'
                self.out.indent()
                self.out.lf() + 'send(&TRUE_BYTE, 1, socketfd, rank);'
                self.out.dedent()
                self.out.lf() + '} else {'
                self.out.indent()
                self.out.lf() + 'send(&FALSE_BYTE, 1, socketfd, rank);'
                self.out.dedent()
                self.out.lf() + "}"
                self.out.dedent();
                self.out.lf() + '}'
                self.out.lf()
                
            else:
                self.out + 'send(' + spec.output_var_name + ', ' + 'header_out[' + spec.counter_name + '] * sizeof(' + spec.type + '), socketfd, rank);'
            
            self.out.dedent().lf() + '}'
            self.out.lf()
            
        self.out.lf().dedent()
        self.out.lf() + '}'
            
        max = self.mapping_from_dtype_to_maximum_number_of_inputvariables.get('string', 0)
        if max != 0:
            self.out.lf() + 'for (int i = 0; i < header_in[HEADER_STRING_COUNT]; i++) {'
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
            self.out.lf() + 'delete[] characters_in;'
            self.out.lf()
            
        max = 0
        dtype = 'string'
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
        if max > 0:
            self.out.lf() + 'delete[] characters_out;'
            self.out.lf()          
            
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
            if dtype == 'string' and max == 0:
                max = 1
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
            self.out.lf()
            self.output_new_statement(
                must_add_type,
                'char *',
                'characters_in',
                max
            )
                
        max = self.mapping_from_dtype_to_maximum_number_of_outputvariables.get(dtype, 0)
        if max == 0:
            max = 1
            
        if max > 0:
            self.out.lf()
            self.output_new_statement(
                must_add_type,
                'char *',
                'characters_out',
                max
            )
                        
        self.out.lf()

            
    def output_new_statement(self, must_add_type, type, var_name, maximum_number_of_inputvariables_of_a_type):
        if must_add_type:
            self.out + type + ' * ' 
            
        self.out + var_name
        self.out + ' = new ' + type
        self.out + '[' + ' max_call_count * ' + maximum_number_of_inputvariables_of_a_type + ']' + ';'
        self.out.lf() 
    
    def output_main(self):
        self.out.lf() + MAIN_FUNCTION_STRING
