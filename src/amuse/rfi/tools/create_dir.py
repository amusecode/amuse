import os

from amuse.support import get_amuse_root_dir
from amuse.support.core import late, print_out
from amuse.support.options import option
from amuse.support.options import OptionalAttributes

interface_file_template = """\
from amuse.community import *

# low level interface class
class {0.name_of_the_community_interface_class}({0.name_of_the_superclass_for_the_community_code_interface_class}):
    
    {0.include_headers_or_modules}
    
    def __init__(self, **keyword_arguments):
        {0.name_of_the_superclass_for_the_community_code_interface_class}.__init__(self, name_of_the_worker="{0.name_of_the_community_code}_worker", **keyword_arguments)
 
# here you must specify the prototypes of the interface functions:
    
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN, unit=None)
        function.addParameter('int_out', dtype='int32', direction=function.OUT, unit=None)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

# optionally, this can be shortened to:
#    @remote_function(can_handle_array=True)
#    def echo_int(int_in='i'):
#        returns (int_out='i')

        
# high level interface class    
class {0.name_of_the_code_interface_class}({0.name_of_the_superclass_for_the_code_interface_class}):

    def __init__(self, **options):
        {0.name_of_the_superclass_for_the_code_interface_class}.__init__(self,  {0.name_of_the_community_interface_class}(**options), **options)

# the following alternative __init__ is appropiate for codes that use an unspecified unit system
# (ie the quantities have dimension but no definite scale) 
#
#    def __init__(self, unit_converter=None, **options):        
#        self.unit_converter=unit_converter
#        {0.name_of_the_superclass_for_the_code_interface_class}.__init__(self,  {0.name_of_the_community_interface_class}(**options), **options)
#
# in this case you also need to use the define_converter below        
        
# typically the high level specification also contains the following:

# the definition of the state model of the code
    def define_state(self, handler):
# for example:
#        handler.set_initial_state('UNINITIALIZED')
#        handler.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
#        handler.add_transition('END', 'STOPPED', 'stop', False)
#        handler.add_transition(
#            'UNINITIALIZED', 'INITIALIZED', 'initialize_code')
#        handler.add_method('STOPPED', 'stop')
        pass

# the definition of any properties
    def define_properties(self, handler):
#        handler.add_property('name_of_the_getter', public_name="name_of_the_property")
        pass

# the definition of the parameters
    def define_parameters(self, handler):
#        handler.add_method_parameter(
#            "name_of_the_getter",
#            "name_of_the_setter",
#            "parameter_name",
#            "description", 
#            default_value = <default value>
#        )
        pass

# the definition of the code data stores, either particle sets:
    def define_particle_sets(self, handler):
#        handler.define_set('particles', 'index_of_the_particle')
#        handler.set_new('particles', 'new_particle')
#        handler.set_delete('particles', 'delete_particle')
#        handler.add_setter('particles', 'set_state')
#        handler.add_getter('particles', 'get_state')
#        handler.add_setter('particles', 'set_mass')
#        handler.add_getter('particles', 'get_mass', names=('mass',))
        pass

# and/or grids:
    def define_grids(self, handler):
#        handler.define_grid('grid',axes_names = ["x", "y"], grid_class=StructuredGrid)
#        handler.set_grid_range('grid', '_grid_range')
#        handler.add_getter('grid', 'get_grid_position', names=["x", "y"])
#        handler.add_getter('grid', 'get_rho', names=["density"])
#        handler.add_setter('grid', 'set_rho', names=["density"])
        pass

# this handles unit conversion if an (optional) unit converter is specified
#    def define_converter(self, handler):
#        if self.unit_converter is not None:
#            handler.set_converter(
#                self.unit_converter.as_converter_from_si_to_generic()
#            )



"""

test_file_template = """\
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
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = {0.name_of_the_interface_code}.o

CODELIB = src/lib{0.name_of_the_community_code}.a

all: {0.name_of_the_community_code}_worker 

clean:
\t$(RM) -rf __pycache__
\t$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
\t$(RM) *~ {0.name_of_the_community_code}_worker worker_code.cc
\tmake -C src clean

distclean: clean
\tmake -C src distclean

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

distclean: clean

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

class CreateADirectoryAndPopulateItWithFiles(OptionalAttributes):
    
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
        return self.amuse_root_dir
        
    @late
    def reference_to_amuse_path(self):
        return os.path.relpath(self.path_of_amuse, self.path_of_the_community_code)
        
    @late 
    def name_of_the_superclass_for_the_community_code_interface_class(self):
        return "CodeInterface"
        
    @late
    def name_of_the_superclass_for_the_code_interface_class(self):
        return "InCodeComponentImplementation"
        
    @late
    def amuse_root_dir(self):
        return get_amuse_root_dir()
        
    @late
    def include_headers_or_modules(self):
        return "include_headers = ['worker_code.h']"
        
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
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = {0.name_of_the_interface_code}.o

CODELIB = src/lib{0.name_of_the_community_code}.a

# needed if code functions are accessed through a module
FCFLAGS+= -I$(realpath ./src)

all: {0.name_of_the_community_code}_worker 

clean:
\t$(RM) -rf __pycache__
\t$(RM) -f *.mod *.so *.o *.pyc worker_code.cc worker_code.h 
\t$(RM) *~ {0.name_of_the_community_code}_worker worker_code.f90
\tmake -C src clean

distclean: clean
\tmake -C src distclean

$(CODELIB):
\tmake -C src all

worker_code.f90: {0.name_of_the_python_module}
\t$(CODE_GENERATOR) --type=f90 interface.py {0.name_of_the_community_interface_class} -o $@

{0.name_of_the_community_code}_worker: worker_code.f90 $(CODELIB) $(OBJS)
\t$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@

%.o: %.f90
\t$(FC) $(FCFLAGS) -c -o $@ $<
"""

code_makefile_template_fortran = """\
MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = lib{0.name_of_the_community_code}.a

CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

clean:
\t$(RM) -f *.o *.a *.mod

distclean: clean

$(CODELIB): $(CODEOBJS)
\t$(RM) -f $@
\t$(AR) $@ $(CODEOBJS)
\t$(RANLIB) $@

%.o: %.f90
\t$(FC) $(FCFLAGS) -c -o $@ $<

"""

code_examplefile_template_fortran = """\
function echo(input)
    integer echo, input
    echo = input
end function
"""

interface_examplefile_template_fortran = """\
module {0.name_of_the_interface_module}


contains

  function echo_int(input, output)
      integer :: echo
      integer :: echo_int
      integer ::  input, output
      output = echo(input)
      echo_int = 0
  end function

end module

"""
class CreateADirectoryAndPopulateItWithFilesForAFortranCode(CreateADirectoryAndPopulateItWithFiles):
        
    @late
    def path_of_the_code_examplefile(self):
        return os.path.join(self.path_of_the_source_code, 'test.f90')
        
    @late
    def path_of_the_interface_examplefile(self):
        return os.path.join(self.path_of_the_community_code, self.name_of_the_interface_code + '.f90')

    @late
    def include_headers_or_modules(self):
        return 'use_modules=["{0}"]'.format(self.name_of_the_interface_module)

    @late
    def name_of_the_interface_module(self):
        return '{0}Interface'.format(self.name_of_the_community_code)

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
            string = interface_examplefile_template_fortran.format(self)
            f.write(string)

