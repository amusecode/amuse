from amuse.support.core import late, print_out

from amuse.legacy import *

import os

interface_file_template = """\
from amuse.legacy import *

class {0.name_of_the_legacy_interface_class}({0.name_of_the_superclass_for_the_legacy_interface_class}):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        {0.name_of_the_superclass_for_the_legacy_interface_class}.__init__(self, **keyword_arguments)
    
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
        
    
class {0.name_of_the_code_interface_class}({0.name_of_the_superclass_for_the_code_interface_class}):

    def __init__(self):
        {0.name_of_the_superclass_for_the_code_interface_class}.__init__(self,  {0.name_of_the_legacy_interface_class}())
    
"""

makefile_template = """\
CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = {0.name_of_the_interface_code}.o

CODELIB = src/lib{0.name_of_the_legacy_code}.a

AMUSE_DIR?={0.reference_to_amuse_path}

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: worker_code 

clean:
\t$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
\t$(RM) *~ worker_code
\tmake -C src clean

$(CODELIB):
\tmake -C src all

worker_code.cc: {0.name_of_the_python_module}
\t$(CODE_GENERATOR) --type=c interface.py {0.name_of_the_legacy_interface_class} -o $@

worker_code.h: {0.name_of_the_python_module}
\t$(CODE_GENERATOR) --type=H interface.py {0.name_of_the_legacy_interface_class} -o $@

worker_code: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
\tmpicxx $(CXXFLAGS) $@.cc $(OBJS) $(CODELIB) -o $@

.cc.o: $<
\t$(CXX) $(CXXFLAGS) -c -o $@ $< 
"""

code_makefile_template = """\
CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = lib{0.name_of_the_legacy_code}.a

CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib

all: $(CODELIB) 


clean:
\t$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
\t$(AR) $@ $(CODEOBJS)
\t$(RANLIB) $@

.cc.o: $<
\t$(CXX) $(CXXFLAGS) -c -o $@ $< 
"""

code_examplefile_template = """\
/*
 * Example function for a code
 */
int echo(int input){
    return input;
}
"""

interface_examplefile_template = """\
extern int echo(int input);

/*
 * Interface code
 */
 
int echo_int(int input, int * output){
    *output = echo(input);
    return 0;
}

"""

class CreateADirectoryAndPopulateItWithFilesForALegacyCode(object):
    
    @late
    def path_of_the_root_directory(self):
        return os.path.dirname(os.path.dirname(__file__))
        
    @late
    def name_of_the_legacy_code(self):
        return self.name_of_the_code_interface_class.lower()
        
    @late
    def name_of_the_python_module(self):
        return 'interface.py'
    
    @late
    def name_of_the_interface_code(self):
        return 'interface'
        
    @late
    def name_of_the_code_interface_class(self):
        return 'MyCode'
        
    @late
    def name_of_the_legacy_interface_class(self):
        return self.name_of_the_code_interface_class + 'Interface'
    
    @late
    def name_of_the_code_directory(self):
        return 'src'
    
    @late
    def path_of_the_legacy_code(self):
        return os.path.join(self.path_of_the_root_directory, self.name_of_the_legacy_code)
        
    @late
    def path_of_the_source_code(self):
        return os.path.join(self.path_of_the_legacy_code, self.name_of_the_code_directory)
        
    @late
    def path_of_the_init_file(self):
        return os.path.join(self.path_of_the_legacy_code, '__init__.py')
        
    @late
    def path_of_the_interface_file(self):
        return os.path.join(self.path_of_the_legacy_code, self.name_of_the_python_module)
        
    @late
    def path_of_the_makefile(self):
        return os.path.join(self.path_of_the_legacy_code, 'Makefile')
    
    @late
    def path_of_the_code_makefile(self):
        return os.path.join(self.path_of_the_source_code, 'Makefile')
        
    @late
    def path_of_the_code_examplefile(self):
        return os.path.join(self.path_of_the_source_code, 'test.cc')
    @late
    def path_of_the_interface_examplefile(self):
        return os.path.join(self.path_of_the_legacy_code, self.name_of_the_interface_code + '.cc')
        
    @late
    def path_of_amuse(self):
        current = os.path.dirname(os.path.dirname(__file__))
        while not os.path.exists(os.path.join(current, 'build.py')):
            current = os.path.dirname(current)
        return current
        
    @late
    def reference_to_amuse_path(self):
        return os.path.relpath(self.path_of_amuse, self.path_of_the_legacy_code)
        
    @late 
    def name_of_the_superclass_for_the_legacy_interface_class(self):
        return LegacyInterface.__name__
        
    @late
    def name_of_the_superclass_for_the_code_interface_class(self):
        return CodeInterface.__name__
        
    def start(self):
        
        self.make_directories()
        self.make_python_files()
        self.make_makefile()
        self.make_example_files()
        
    def make_directories(self):
        os.mkdir(self.path_of_the_legacy_code)
        os.mkdir(self.path_of_the_source_code)
        
    def make_python_files(self):
        with open(self.path_of_the_init_file, "w") as f:
            f.write("# generated file")
        
        with open(self.path_of_the_interface_file, "w") as f:
            string = interface_file_template.format(self)
            f.write(string)
        
    def make_makefile(self):
        
        with open(self.path_of_the_makefile, "w") as f:
            string = makefile_template.format(self)
            f.write(string)
            
    def make_example_files(self):
        with open(self.path_of_the_code_makefile, "w") as f:
            string = code_makefile_template.format(self)
            f.write(string)
            
        with open(self.path_of_the_code_examplefile, "w") as f:
            string = code_examplefile_template
            f.write(string)
        
        with open(self.path_of_the_interface_examplefile, "w") as f:
            string = interface_examplefile_template
            f.write(string)
        
        
