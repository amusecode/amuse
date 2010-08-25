from amuse.support.legacy.core import *

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.legacy import channel

from amuse.test.amusetest import TestWithMPI
from amuse.support.legacy import create_c
#cello
from amuse.support.legacy import create_fortran

from amuse.support.legacy import stopping_conditions
from amuse.support.interface import CodeInterface

import subprocess
import os

codestring = """
#include <stopcond.h>
#ifdef __cplusplus
extern "C" {
#endif
int initialize_code()
{
    
    // AMUSE STOPPING CONDITIONS SUPPORT
    supported_conditions = COLLISION_DETECTION_BITMAP | PAIR_DETECTION_BITMAP | TIMEOUT_DETECTION_BITMAP;
    // -----------------------
    return 0;
}
#ifdef __cplusplus
}
#endif
"""

codestringF = """
      FUNCTION initialize_code()
      IMPLICIT NONE
      include "../../lib/stopcond/stopcond.inc"
      INTEGER :: initialize_code
      INTEGER :: set_supported_conditions
      INTEGER :: return
      initialize_code = 0
      return = set_supported_conditions(IOR(COLLISION_DETECTION_BITMAP, PAIR_DETECTION_BITMAP))
      RETURN
      END

"""

class ForTestingInterface(LegacyInterface, stopping_conditions.StoppingConditionInterface):
    def __init__(self, exefile, **options):
        LegacyInterface.__init__(self, exefile, **options)

    include_headers = ['stopcond.h']

    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function     

    @legacy_function
    def reset_stopping_conditions():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        function.can_handle_array = False
        return function     
        
    @legacy_function
    def next_index_for_stopping_condition():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.can_handle_array = False
        return function  
        
    @legacy_function
    def set_stopping_condition_info():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_condition', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
    
    @legacy_function
    def set_stopping_condition_particle_index():
        function = LegacyFunctionSpecification()  
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_condition', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  

class ForTesting(CodeInterface):
    def __init__(self, exefile, **options):
        self.stopping_conditions = stopping_conditions.StoppingConditions(self)
        CodeInterface.__init__(self, ForTestingInterface(exefile, **options), **options)
    
    def define_methods(self, object):
        self.stopping_conditions.define_methods(object)

        
class TestInterface(TestWithMPI):
    def cxx_compile(self, objectname, string):
  
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.cc'
        if os.path.exists(objectname):
            os.remove(objectname)
        with open(sourcename, "w") as f:
            f.write(string)
        
        rootdir = self.get_amuse_root_dir()
        arguments = ["mpicxx", "-I",rootdir + "/lib/stopcond", "-c",  "-o", objectname, sourcename]
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if not os.path.exists(objectname): # or process.poll() == 1:
            raise Exception("Could not compile {0}, error = {1} ({2})".format(objectname, stderr, ' '.join(arguments)))
    
    def c_build(self, exename, objectnames):
        rootdir = self.get_amuse_root_dir()
        
        arguments = ["mpicxx"]
        arguments.extend(objectnames)
        arguments.append("-o")
        arguments.append(exename)
        arguments.extend(["-L"+rootdir+"/lib/stopcond","-lstopcond"])
        print ' '.join(arguments)
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))
    
    def build_worker(self):
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"code.o")
        interfacefile = os.path.join(path,"interface.o")
        self.exefile = os.path.join(path,"c_worker")
        
        self.cxx_compile(codefile, codestring)
        
        uc = create_c.MakeACHeaderStringOfAClassWithLegacyFunctions()
        uc.class_with_legacy_functions = ForTestingInterface
        uc.make_extern_c = True
        uc.ignore_functions_from = [stopping_conditions.StoppingConditionInterface]
        header =  uc.result

        uc = create_c.MakeACStringOfAClassWithLegacyFunctions()

        uc.class_with_legacy_functions = ForTestingInterface
        code =  uc.result

        string = '\n\n'.join([header, code])
        
        #print string
        
        self.cxx_compile(interfacefile, string)
        self.c_build(self.exefile, [interfacefile, codefile] )
    
    def setUp(self):
        super(TestInterface, self).setUp()
        print "building"
        self.build_worker()
        
    def test1(self):
        print self.exefile
        instance = ForTestingInterface(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        next = instance.next_index_for_stopping_condition()
        instance.stop()
        self.assertEquals(next, 1)
        
    def test2(self):
        instance = ForTesting(self.exefile, redirection = "none") #, debugger = "xterm")
        instance.initialize_code()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_supported())
        self.assertTrue(instance.stopping_conditions.collision_detection.is_supported())
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        instance.stop()
        
    def test3(self):
        instance = ForTesting(self.exefile, redirection = "none") #, debugger = "xterm")
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.enable()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.disable()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stop()
        
    def test4(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())

        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
         
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.stop()
        
    def test5(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())

        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, 11)
        instance.set_stopping_condition_particle_index(next, 1, 12)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEquals(12, instance.get_stopping_condition_particle_index(next, 1))
        instance.stop()

class TestInterfaceF(TestWithMPI):
    def f90_compile(self, objectname, string):
        root, ext = os.path.splitext(objectname)
        sourcename = root + '.f90'
        if os.path.exists(objectname):
            os.remove(objectname)
        with open(sourcename, "w") as f:
            f.write(string)
        
        rootdir = self.get_amuse_root_dir()
        arguments = ["mpif90", "-I",rootdir + "/lib/stopcond", "-c",  "-o", objectname, sourcename]
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if not os.path.exists(objectname): # or process.poll() == 1:
            raise Exception("Could not compile {0}, error = {1} ({2})".format(objectname, stderr, ' '.join(arguments)))
    
    def f90_build(self, exename, objectnames):
        rootdir = self.get_amuse_root_dir()
        
        arguments = ["mpif90"]
        arguments.extend(objectnames)
        arguments.append("-o")
        arguments.append(exename)
        arguments.extend(["-L"+rootdir+"/lib/stopcond","-lstopcond"])
        print 'build command:'
        print ' '.join(arguments)
        process = subprocess.Popen(
            arguments,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise Exception("Could not build {0}, error = {1}".format(exename, stderr))
    
    def build_worker(self):
        path = os.path.abspath(self.get_path_to_results())
        codefile = os.path.join(path,"codef90.o")
        interfacefile = os.path.join(path,"interfacef90.o")
        self.exefile = os.path.join(path,"fortran_worker")
        
        self.f90_compile(codefile, codestringF)
        
        uf = create_fortran.MakeAFortranStringOfAClassWithLegacyFunctions()
        uf.class_with_legacy_functions = ForTestingInterface
        string =  uf.result
        
        #print string
        
        self.f90_compile(interfacefile, string)
        self.f90_build(self.exefile, [interfacefile, codefile] )
    
    def setUp(self):
        super(TestInterfaceF, self).setUp()
        print "building"
        self.build_worker()

    def test1(self):
        print self.exefile
        instance = ForTestingInterface(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        next = instance.next_index_for_stopping_condition()
        instance.stop()
        self.assertEquals(next, 1)

    def test2(self):
        instance = ForTesting(self.exefile, redirection = "none") #, debugger = "xterm")
        instance.initialize_code()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_supported())
        self.assertTrue(instance.stopping_conditions.collision_detection.is_supported())
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        instance.stop()

    def test3(self):
        instance = ForTesting(self.exefile, redirection = "none") #, debugger = "xterm")
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.enable()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.disable()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stop()

                
    def test4(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()

        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
    
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
         
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.stop()
    
    def test5(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        
        instance.set_stopping_condition_info(next,instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, 11)
        instance.set_stopping_condition_particle_index(next, 1, 12)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEquals(12, instance.get_stopping_condition_particle_index(next, 1))
        instance.stop()
    
