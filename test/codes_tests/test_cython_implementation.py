from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.test.amusetest import TestCase
from amuse.support import exceptions
from amuse.support import options

import subprocess
import os
import time

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
        

class TestCythonImplementationInterface(test_c_implementation.TestCImplementationInterface):

    def setUp(self):
        self.skip_if_no_cython()
        super(TestCythonImplementationInterface, self).setUp()


    def tearDown(self):
        pass

    def test22(self):
        self.skip("this test uses mpi internals, skip here")
                     
    def skip_if_no_cython(self):
        
        if not hasattr(config.compilers, 'cython') or len(config.compilers.cython) == 0:
            self.skip("cython not configured")
            
        process = subprocess.Popen(
            [config.compilers.cython,'-V'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        output, error = process.communicate()
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
        
        self.c_compile(codefile, test_c_implementation.codestring)
        
        
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
        
        process = subprocess.Popen(
           # [config.compilers.cython, '--embed', sourcename, '-o', cname],
            [config.compilers.cython,  sourcename, '-o', cname],
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        
        if process.returncode == 0:
            self.wait_for_file(cname)
        
        if process.returncode != 0 or not os.path.exists(cname):
            print "Could not cythonize {0}, error = {1}".format(sourcename, stderr)
            raise Exception("Could not cythonize {0}, error = {1}".format(sourcename, stderr))
        
        with open(cname, "r") as f:
            string = f.read()
            
        self.c_pythondev_compile(interfacefile, string)
        # self.c_pythondev_build(self.exefile, [interfacefile, codefile] )
        self.c_pythondev_buildso(self.sofile,  [interfacefile, codefile] )



    def c_pythondev_compile(self, objectname, string):
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.c'
        
        if os.path.exists(objectname):
            os.remove(objectname)
        
        with open(sourcename, "w") as f:
            f.write(string)
            
        mpicc = self.get_mpicc_name()
        arguments = [mpicc]
        arguments.extend(self.get_mpicc_flags().split())
        arguments.extend(["-I", "lib/stopcond","-I", "lib/amuse_mpi",  "-fPIC", "-c",  "-o", objectname, sourcename])
        arguments.extend(shlex.split(config.compilers.pythondev_cflags))
            
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(objectname)
        
        if process.returncode != 0 or not os.path.exists(objectname):
            print "Could not compile {0}, error = {1}".format(objectname, stderr)
            raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))


    def c_pythondev_build(self, exename, objectnames):
        
        if os.path.exists(exename):
            os.remove(exename)
            
        mpicxx = self.get_mpicxx_name()
        arguments = [mpicxx]
        arguments.extend(objectnames)
        arguments.extend(shlex.split(config.compilers.pythondev_ldflags))
        arguments.append("-o")
        arguments.append(exename)
        
        if 'LIBS' in os.environ:
            libs = os.environ['LIBS'].split()
            arguments.extend(libs)
            
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(exename)
            
        if process.returncode != 0 or not (os.path.exists(exename) or os.path.exists(exename+'.exe')):
            print "Could not compile {0}, error = {1}".format(exename, stderr)
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))

        print stdout
        print stderr


    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        
    def test14(self):
        self.skip("needs python support")
                     
    def c_pythondev_buildso(self, soname, objectnames):
        
        if os.path.exists(soname):
            os.remove(soname)
            
        mpicxx = self.get_mpicxx_name()
        arguments = [mpicxx]
        arguments.extend(objectnames)
        arguments.extend(shlex.split(config.compilers.pythondev_ldflags))
        arguments.append("-shared")
        arguments.append("-o")
        arguments.append(soname)
        arguments.append("-Llib/amuse_mpi")
        arguments.append("-lamuse_mpi")
        
        if 'LIBS' in os.environ:
            libs = os.environ['LIBS'].split()
            arguments.extend(libs)
            
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(soname)
            
        if process.returncode != 0 or not (os.path.exists(soname)):
            print "Could not compile {0}, error = {1}".format(soname, stderr)
            raise Exception("Could not build {0}, error = {1}".format(soname, stderr))

        print stdout
        print stderr




    def c_compile(self, objectname, string):
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.c'
        
        if os.path.exists(objectname):
            os.remove(objectname)
        
        with open(sourcename, "w") as f:
            f.write(string)
            
        mpicc = self.get_mpicc_name()
        arguments = [mpicc]
        arguments.extend(self.get_mpicc_flags().split())
        arguments.extend(["-I", "lib/stopcond", "-fPIC", "-c",  "-o", objectname, sourcename])
            
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            self.wait_for_file(objectname)
        
        if process.returncode != 0 or not os.path.exists(objectname):
            print "Could not compile {0}, error = {1}".format(objectname, stderr)
            raise Exception("Could not compile {0}, error = {1}".format(objectname, stderr))


