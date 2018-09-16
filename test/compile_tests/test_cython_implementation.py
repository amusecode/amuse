from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.test.amusetest import TestCase
from amuse.test import compile_tools
from amuse.support import exceptions
from amuse.support import options

import os
import time
import sys


from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi.tools import create_c
from amuse.rfi import channel
from amuse.rfi.core import *

import test_c_implementation
from amuse.rfi.tools import create_cython

codestring = """
#include <stdio.h>
#include <stdlib.h>

int echo_int(int int_in, int * int_out) {
    *int_out = int_in;
    if(int_in < 0) {
        return -1;
    } else {
        return 0;
    }
}"""


import shlex
from amuse.test.amusetest import get_amuse_root_dir

import test_fortran_implementation

from amuse.rfi.tools import create_fortran

class ForTestingInterface(CodeInterface):
    
    def __init__(self, exefile, **options):
        CodeInterface.__init__(self, exefile, **options)

    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function 
            
    
class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)
    
    def define_methods(self, object):
        object.add_method(
            'echo_int',
            (units.m,)
            ,
            (
                units.m,
                object.ERROR_CODE,
            )
        )

            
class TestCreateCython(TestCase):
    
    def test1(self):
        uc = create_cython.GenerateACythonSourcecodeStringFromASpecificationClass()
        uc.specification_class = test_c_implementation.ForTestingInterface
        uc.needs_mpi = True
        code =  uc.result
        #print code
        self.assertTrue('cdef extern from "worker_code.h":' in code)
        self.assertTrue('int c_echo_int "echo_int" (int, int *)' in code)
        

    def test2(self):
        uc = create_cython.GenerateAFortranInterfaceStringOfAFunctionSpecification()
        uc.specification = test_c_implementation.ForTestingInterface.echo_int.specification
        uc.needs_mpi = True
        uc.start()
        code =  uc.result
        print '<c>' + code + '</c>'
        self.assertTrue('INTEGER(kind = c_int) function c_3f33a9ce(int_in, &\n &int_out) &\n & result(rrreeesss) &\n & bind(c, name = "ci_echo_int")' in code)
        self.assertTrue('  INTEGER(kind = c_int), intent(in), value :: int_in' in code)
        self.assertTrue('  INTEGER(kind = c_int), intent(out) :: int_out' in code)
        self.assertTrue('  INTEGER :: echo_int' in code)
        self.assertTrue('  rrreeesss = echo_int(int_in, &\n &int_out)' in code)
        self.assertTrue('END FUNCTION c_3f33a9ce' in code)
        



    def test3(self):
        uc = create_cython.GenerateAFortranInterfaceStringOfAFunctionSpecification()
        uc.specification = test_c_implementation.ForTestingInterface.echo_string.specification
        uc.needs_mpi = True
        uc.start()
        code =  uc.result
        print '<c>' + code + '</c>'
        self.assertTrue('INTEGER(kind = c_int) function c_64452f1a(string_in, &\n &string_out) &\n & result(rrreeesss) &\n & bind(c, name = "ci_echo_string")' in code)
        self.assertTrue('  type(C_ptr), intent(in), value :: string_in' in code)
        self.assertTrue('  character(len=4096) :: string_string_in' in code)
        self.assertTrue('  call C_F_string_ptr(string_in, string_string_in)' in code)
        self.assertTrue('  call C_F_string_ptr(string_in, string_string_in)' in code)
        self.assertTrue('  rrreeesss = echo_string(string_string_in, &\n &string_string_out)' in code)
        self.assertTrue('     string_buffer1 = C_malloc(sz)' in code)
        self.assertTrue('  call F_C_string_ptr(trim(string_string_out), string_buffer1, 4096)' in code)
        self.assertTrue('END FUNCTION c_64452f1a' in code)
        




class TestCythonImplementationInterface(test_c_implementation.TestCImplementationInterface):

    def setUp(self):
        self.skip_if_no_cython()
        super(test_c_implementation.TestCImplementationInterface, self).setUp()
        print "building...",
        self.check_can_compile_modules()
        try:
            self.build_worker()
        except Exception as ex:
            print ex
            raise
        print "done"


    def tearDown(self):
        pass

    def test22(self):
        self.skip("this test uses mpi internals, skip here")
                     
    def skip_if_no_cython(self):
        
        if sys.hexversion > 0x03000000:
            self.skip("no cython for python 3.0")

        if not hasattr(config.compilers, 'cython') or len(config.compilers.cython) == 0:
            self.skip("cython not configured")
            
        process, output, error = compile_tools.open_subprocess([config.compilers.cython, '-V'])
        if process.returncode:
            self.skip("cython not available")
        
        error = error.strip()[0:7].decode('utf-8')
        if not error.startswith('Cython'):
            self.skip("cython not available")
        
            
        
        
            

    def build_worker(self):
        
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        interfacefile = os.path.join(path,"interface.o")
        headerfile = os.path.join(path,"worker_code.h")
        self.sofile = os.path.join(path,"interface.so")
        self.exefile = os.path.join(path,"c_worker")
        
        compile_tools.c_compile(codefile, test_c_implementation.codestring,
                                extra_args=["-fPIC"])
        
        uc = create_c.GenerateACHeaderStringFromASpecificationClass()
        uc.specification_class = test_c_implementation.ForTestingInterface
        uc.needs_mpi = False
        header =  uc.result
        
        
        with open(headerfile, "w") as f:
            f.write(header)
        
        
        root, ext = os.path.splitext(interfacefile)
        sourcename = root + '.pyx'
        cname = root + '.c'
        
        uc = create_cython.GenerateACythonSourcecodeStringFromASpecificationClass()
        uc.specification_class = test_c_implementation.ForTestingInterface
        uc.needs_mpi = True
        code =  uc.result

        


        
        with open(sourcename, "w") as f:
            f.write(code)


        uc = create_cython.GenerateACythonStartScriptStringFromASpecificationClass()
        uc.specification_class = test_c_implementation.ForTestingInterface
        uc.needs_mpi = True
        script =  uc.result
        
        with open(self.exefile, "w") as f:
            f.write(script)

        os.chmod(self.exefile, 0777)
        
        import mpi4py
        process, stdout, stderr = compile_tools.open_subprocess([config.compilers.cython, 
        '-I',
        mpi4py.get_include(),
         sourcename, '-o', cname])

        
        if process.returncode == 0:
            compile_tools.wait_for_file(cname)
        
        if process.returncode != 0 or not os.path.exists(cname):
            print "Could not cythonize {0}, error = {1}".format(sourcename, stderr)
            raise Exception("Could not cythonize {0}, error = {1}".format(sourcename, stderr))
        
        with open(cname, "r") as f:
            string = f.read()
            
        compile_tools.c_pythondev_compile(interfacefile, string)
        # self.c_pythondev_build(self.exefile, [interfacefile, codefile] )
        compile_tools.c_pythondev_buildso(self.sofile,  [interfacefile, codefile] )

    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        
    def test14(self):
        self.skip("needs python support")


class TestCythonFortranImplementationInterface(test_fortran_implementation.TestInterface):

    def setUp(self):
        self.skip_if_no_cython()
        super(test_fortran_implementation.TestInterface, self).setUp()
        print "building"
        self.check_can_compile_modules()
        self.build_worker()


    def tearDown(self):
        pass

    def test22(self):
        self.skip("this test uses mpi internals, skip here")
                     
    def skip_if_no_cython(self):

        if sys.hexversion > 0x03000000:
            self.skip("no cython for python 3.0")
        
        if not hasattr(config.compilers, 'cython') or len(config.compilers.cython) == 0:
            self.skip("cython not configured")
            
        process, output, error = compile_tools.open_subprocess([config.compilers.cython, '-V'])
        if process.returncode:
            self.skip("cython not available")
        
        error = error.strip()[0:7].decode('utf-8')
        if not error.startswith('Cython'):
            self.skip("cython not available")

    def build_worker(self):
        
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        interfacefile = os.path.join(path,"interface.o")
        headerfile = os.path.join(path,"worker_code.h")
        self.sofile = os.path.join(path,"interface.so")
        self.interfacec_o_file = os.path.join(path,"interfacec.o")
        self.exefile = os.path.join(path,"c_worker")
        
        compile_tools.fortran_compile(codefile, test_fortran_implementation.codestring,
                                      extra_args=["-fPIC"])
        
        
        uc = create_c.GenerateACHeaderStringFromASpecificationClass()
        uc.specification_class = test_fortran_implementation.ForTestingInterface
        uc.needs_mpi = False
        header =  uc.result
        
        
        with open(headerfile, "w") as f:
            f.write(header)
        
        
        root, ext = os.path.splitext(interfacefile)
        sourcename = root + '.pyx'
        cname = root + '.c'
        
        uc = create_cython.GenerateACythonSourcecodeStringFromASpecificationClass()
        uc.specification_class = test_fortran_implementation.ForTestingInterface
        uc.function_name_prefix = "ci_"
        uc.needs_mpi = True
        code =  uc.result
        
        with open(sourcename, "w") as f:
            f.write(code)


        uc = create_cython.GenerateACythonStartScriptStringFromASpecificationClass()
        uc.specification_class = test_fortran_implementation.ForTestingInterface
        uc.needs_mpi = True
        script =  uc.result
        
        with open(self.exefile, "w") as f:
            f.write(script)

        os.chmod(self.exefile, 0777)

        import mpi4py
        process, stdout, stderr = compile_tools.open_subprocess([config.compilers.cython, 
        '-I',
        mpi4py.get_include(),
         sourcename, '-o', cname])

        if process.returncode == 0:
            compile_tools.wait_for_file(cname)
        
        if process.returncode != 0 or not os.path.exists(cname):
            print "Could not cythonize {0}, error = {1}".format(sourcename, stderr)
            raise Exception("Could not cythonize {0}, error = {1}".format(sourcename, stderr))
        
        with open(cname, "r") as f:
            string = f.read()


        
        
        uc = create_cython.GenerateAFortranInterfaceSourcecodeStringFromASpecificationClass()
        uc.specification_class = test_fortran_implementation.ForTestingInterface
        uc.function_name_prefix = "ci_"
        uc.needs_mpi = False
        code =  uc.result

        compile_tools.fortran_compile(self.interfacec_o_file, code,
                                      extra_args=["-fPIC"])

        compile_tools.c_pythondev_compile(interfacefile, string)
        compile_tools.fortran_pythondev_buildso(self.sofile,  [interfacefile, codefile, self.interfacec_o_file] )

    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
