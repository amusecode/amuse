from amuse.support.core import late, print_out

import os

interface_file_template = """\
from amuse.community import *

class {0.name_of_the_community_interface_class}({0.name_of_the_superclass_for_the_community_code_interface_class}):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        {0.name_of_the_superclass_for_the_community_code_interface_class}.__init__(self, name_of_the_worker="{0.name_of_the_community_code}_worker", **keyword_arguments)
    
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
        {0.name_of_the_superclass_for_the_code_interface_class}.__init__(self,  {0.name_of_the_community_interface_class}())
    
"""

test_file_template = """\
from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from {0.name_for_import_of_the_interface_module} import {0.name_of_the_community_interface_class}
from {0.name_for_import_of_the_interface_module} import {0.name_of_the_code_interface_class}

class {0.name_of_the_community_interface_class}Tests(TestWithMPI):
    
    def test1(self):
        instance = {0.name_of_the_community_interface_class}()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
"""

makefile_template_cxx = """\
# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?={0.reference_to_amuse_path}
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = {0.name_of_the_interface_code}.o

CODELIB = src/lib{0.name_of_the_community_code}.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: {0.name_of_the_community_code}_worker 

clean:
\t$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
\t$(RM) *~ {0.name_of_the_community_code}_worker worker_code.cc
\tmake -C src clean

$(CODELIB):
\tmake -C src all

worker_code.cc: {0.name_of_the_python_module}
\t$(CODE_GENERATOR) --type=c interface.py {0.name_of_the_community_interface_class} -o $@

worker_code.h: {0.name_of_the_python_module}
\t$(CODE_GENERATOR) --type=H interface.py {0.name_of_the_community_interface_class} -o $@

{0.name_of_the_community_code}_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
\t$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
\t$(CXX) $(CXXFLAGS) -c -o $@ $< 
"""

code_makefile_template_cxx = """\
CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = lib{0.name_of_the_community_code}.a

CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 


clean:
\t$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
\t$(RM) -f $@
\t$(AR) $@ $(CODEOBJS)
\t$(RANLIB) $@

.cc.o: $<
\t$(CXX) $(CXXFLAGS) -c -o $@ $< 
"""

code_examplefile_template_cxx = """\
/*
 * Example function for a code
 */
int echo(int input){
    return input;
}
"""

interface_examplefile_template_cxx = """\
extern int echo(int input);

/*
 * Interface code
 */
 
int echo_int(int input, int * output){
    *output = echo(input);
    return 0;
}

"""

class CreateADirectoryAndPopulateItWithFiles(object):
    
    @late
    def path_of_the_root_directory(self):
        return os.path.dirname(os.path.dirname(__file__))
        
    @late
    def name_of_the_community_code(self):
        return self.name_of_the_code_interface_class.lower()
        
    @late
    def name_of_the_python_module(self):
        return 'interface.py'
        
    @late
    def name_of_the_test_module(self):
        return 'test_{0}.py'.format(self.name_of_the_community_code)
    
    @late
    def name_of_the_interface_code(self):
        return 'interface'
        
    @late
    def name_of_the_code_interface_class(self):
        return 'MyCode'
        
    @late
    def name_of_the_community_interface_class(self):
        return self.name_of_the_code_interface_class + 'Interface'
    
    @late
    def name_of_the_code_directory(self):
        return 'src'
    
    @late
    def name_for_import_of_the_interface_module(self):
        return '.' + self.name_of_the_python_module[:-3]
        
    @late
    def path_of_the_community_code(self):
        return os.path.join(self.path_of_the_root_directory, self.name_of_the_community_code)
        
    @late
    def path_of_the_source_code(self):
        return os.path.join(self.path_of_the_community_code, self.name_of_the_code_directory)
        
    @late
    def path_of_the_init_file(self):
        return os.path.join(self.path_of_the_community_code, '__init__.py')
        
    @late
    def path_of_the_interface_file(self):
        return os.path.join(self.path_of_the_community_code, self.name_of_the_python_module)
    
    @late
    def path_of_the_test_file(self):
        return os.path.join(self.path_of_the_community_code, self.name_of_the_test_module)
        
    @late
    def path_of_the_makefile(self):
        return os.path.join(self.path_of_the_community_code, 'Makefile')
    
    @late
    def path_of_the_code_makefile(self):
        return os.path.join(self.path_of_the_source_code, 'Makefile')
        
    @late
    def path_of_the_code_examplefile(self):
        raise NotImplementedError()
        
    @late
    def path_of_the_interface_examplefile(self):
        raise NotImplementedError()
        
    @late
    def path_of_amuse(self):
        current = os.path.dirname(os.path.dirname(__file__))
        while not os.path.exists(os.path.join(current, 'build.py')):
            current = os.path.dirname(current)
        return current
        
    @late
    def reference_to_amuse_path(self):
        return os.path.relpath(self.path_of_amuse, self.path_of_the_community_code)
        
    @late 
    def name_of_the_superclass_for_the_community_code_interface_class(self):
        return "CodeInterface"
        
    @late
    def name_of_the_superclass_for_the_code_interface_class(self):
        return "InCodeComponentImplementation"
        
    def start(self):
        
        self.make_directories()
        self.make_python_files()
        self.make_makefile()
        self.make_example_files()
        
        
    def make_directories(self):
        os.mkdir(self.path_of_the_community_code)
        os.mkdir(self.path_of_the_source_code)
        
    def make_python_files(self):
        with open(self.path_of_the_init_file, "w") as f:
            f.write("# generated file")
        
        with open(self.path_of_the_interface_file, "w") as f:
            string = interface_file_template.format(self)
            f.write(string)
        
        with open(self.path_of_the_test_file, "w") as f:
            string = test_file_template.format(self)
            f.write(string)
        
    def make_makefile(self):
        pass
            
    def make_example_files(self):
        pass
        
        
class CreateADirectoryAndPopulateItWithFilesForACCode(CreateADirectoryAndPopulateItWithFiles):
   
    @late
    def path_of_the_code_examplefile(self):
        return os.path.join(self.path_of_the_source_code, 'test.cc')
        
    @late
    def path_of_the_interface_examplefile(self):
        return os.path.join(self.path_of_the_community_code, self.name_of_the_interface_code + '.cc')
            
    def make_makefile(self):
        
        with open(self.path_of_the_makefile, "w") as f:
            string = makefile_template_cxx.format(self)
            f.write(string)
            
    def make_example_files(self):
        with open(self.path_of_the_code_makefile, "w") as f:
            string = code_makefile_template_cxx.format(self)
            f.write(string)
            
        with open(self.path_of_the_code_examplefile, "w") as f:
            string = code_examplefile_template_cxx
            f.write(string)
        
        with open(self.path_of_the_interface_examplefile, "w") as f:
            string = interface_examplefile_template_cxx
            f.write(string)



makefile_template_fortran = """\
# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?={0.reference_to_amuse_path}
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

FFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = {0.name_of_the_interface_code}.o

CODELIB = src/lib{0.name_of_the_community_code}.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: {0.name_of_the_community_code}_worker 

clean:
\t$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
\t$(RM) *~ worker_code worker_code.f90
\tmake -C src clean

$(CODELIB):
\tmake -C src all

worker_code.f90: {0.name_of_the_python_module}
\t$(CODE_GENERATOR) --type=f90 interface.py {0.name_of_the_community_interface_class} -o $@

{0.name_of_the_community_code}_worker: worker_code.f90 $(CODELIB) $(OBJS)
\t$(MPIFC) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

%.o: %.f90
\t$(FC) $(FFLAGS) -c -o $@ $<
"""

code_makefile_template_fortran = """\
MPIFC ?= mpif90
FC      = $(MPIFC)

FFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = lib{0.name_of_the_community_code}.a

CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

clean:
\t$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
\t$(RM) -f $@
\t$(AR) $@ $(CODEOBJS)
\t$(RANLIB) $@

%.o: %.f90
\t$(FC) $(FFLAGS) -c -o $@ $<

"""

code_examplefile_template_fortran = """\
FUNCTION echo(input)
    INTEGER echo, input
    echo = input
END FUNCTION
"""

interface_examplefile_template_fortran = """\
FUNCTION echo_int(input, output)
    INTEGER echo
    INTEGER echo_int
    INTEGER input, output
    output = echo(input)
    echo_int = 0
END FUNCTION

"""
class CreateADirectoryAndPopulateItWithFilesForAFortranCode(CreateADirectoryAndPopulateItWithFiles):
   
    @late
    def path_of_the_code_examplefile(self):
        return os.path.join(self.path_of_the_source_code, 'test.f90')
        
    @late
    def path_of_the_interface_examplefile(self):
        return os.path.join(self.path_of_the_community_code, self.name_of_the_interface_code + '.f90')
            
    def make_makefile(self):
        
        with open(self.path_of_the_makefile, "w") as f:
            string = makefile_template_fortran.format(self)
            f.write(string)
            
    def make_example_files(self):
        with open(self.path_of_the_code_makefile, "w") as f:
            string = code_makefile_template_fortran.format(self)
            f.write(string)
            
        with open(self.path_of_the_code_examplefile, "w") as f:
            string = code_examplefile_template_fortran
            f.write(string)
        
        with open(self.path_of_the_interface_examplefile, "w") as f:
            string = interface_examplefile_template_fortran
            f.write(string)

