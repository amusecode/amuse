from amuse.legacy.support.core import *

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.legacy.support import python_code
from amuse.legacy.support import channel

from legacy_support import TestWithMPI

import parser

class ForTestingInterface(LegacyPythonInterface):
    
    def __init__(self):
        LegacyPythonInterface.__init__(self, implementation_factory = ForTestingImplementation)
        
    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True 
        function.id = 10
        return function
        
    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 11
        return function    
        
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 12
        return function     
        
    @legacy_function
    def echo_double():
        function = LegacyFunctionSpecification()  
        function.addParameter('double_in', dtype='float64', direction=function.IN)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 13
        return function  
          
    @legacy_function
    def echo_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.addParameter('string_out', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 14
        return function  
          
    @legacy_function
    def echo_strings():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_inout1', dtype='string', direction=function.INOUT)
        function.addParameter('string_inout2', dtype='string', direction=function.INOUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        function.id = 15
        return function    

    
    
    
    
    
class ForTestingImplementation(object):
    
    def __init__(self):
        self.masses = [0.0] * 100
        
    def get_mass(self, index_of_the_particle,  mass):
        try:
            mass.value = self.masses[index_of_the_particle]
            return 0
        except:
            return -1
        
    def set_mass(self, index_of_the_particle,  mass):
        try:
            self.masses[index_of_the_particle] = mass
            return 0
        except:
            return -1
            
    def echo_int(self, int_in, int_out):
        int_out.value = int_in
        return 0
            
    def echo_double(self, double_in, double_out):
        double_out.value = double_in
        return 0
        
    def echo_string(self, string_in, string_out):
        string_out.value = string_in
        return 0
        
    def echo_strings(self, string_inout1, string_inout2):
        string_inout1.value = string_inout1.value[::-1]
        string_inout2.value = string_inout2.value[::-1]
        return 0
       

class TestInterface(TestWithMPI):
    
    def test1(self):
        script_string = ForTestingInterface.new_executable_script_string_for(ForTestingImplementation)
        
        self.assertTrue(script_string.find('syspath = (') > 0)
        self.assertTrue(script_string.find('ForTestingInterface') > 0)
        self.assertTrue(script_string.find('ForTestingImplementation') > 0)
        self.assertTrue(script_string.find('test_python_implementation') > 0)
        self.assertTrue(script_string.find('PythonImplementation(instance, ForTestingInterface)')>0)
        try:
            st = compile(script_string, 'test.py', 'exec')
        except SyntaxError, ex:
            self.fail("Compilation error {0}".format(ex))
            
    def test2(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.Message(10, 1)
        input_message.ints = [1]
        
        output_message = python_code.Message(10, 1)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 1)
        self.assertEquals(len(output_message.doubles), 1)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(output_message.doubles[0], 0.0)
        
    def test3(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.Message(11, 1)
        input_message.ints = [1]
        input_message.doubles = [12.0]
        
        output_message = python_code.Message(10, 1)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 1)
        self.assertEquals(len(output_message.doubles), 0)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(implementation.masses[1], 12.0)
        
    
    def test4(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.Message(11, 4)
        input_message.ints = [1,2,3,4]
        input_message.doubles = [12.0,13.0,14.0,15.0]
        
        output_message = python_code.Message(10, 4)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 4)
        self.assertEquals(len(output_message.doubles), 0)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(output_message.ints[3], 0)
        self.assertEquals(implementation.masses[1], 12.0)
        self.assertEquals(implementation.masses[2], 13.0)
        self.assertEquals(implementation.masses[3], 14.0)
        self.assertEquals(implementation.masses[4], 15.0)
        
    def test5(self):
        x = ForTestingInterface()
        
        error = x.set_mass(1, 10.0)
        self.assertEquals(error, 0)
        
        answer, error = x.get_mass(1)
        self.assertEquals(error, 0)
        self.assertEquals(answer, 10.0)
        
        
        del x
        
        
    def test6(self):
        x = ForTestingInterface()
        
        errors = x.set_mass([1,2], [10.0,11.0])
        self.assertEquals(errors[0], 0)
        self.assertEquals(errors[1], 0)
        
        answer, errors = x.get_mass([1,2])
        self.assertEquals(errors[0], 0)
        self.assertEquals(answer[0], 10.0)
        self.assertEquals(answer[1], 11.0)
        
        del x
        
    def test7(self):
        x = ForTestingInterface()
        
        int_out, error = x.echo_int(20)
        self.assertEquals(error, 0)
        self.assertEquals(int_out, 20)
        
        del x
        
    
    def test8(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.Message(12, 1)
        input_message.ints = [20]
        
        output_message = python_code.Message(10, 1)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 2)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(output_message.ints[1], 20)
        
    
    def test9(self):
        x = ForTestingInterface()
        channel.MessageChannel.DEBUGGER = None
        string_out, error = x.echo_string("1234567")
        self.assertEquals(error, 0)
        self.assertEquals(string_out[0], "1234567")
        
        del x
        
    def test10(self):
        x = ForTestingInterface()
        string_out, error = x.echo_string(["aaaaa", "bbbb"])
        self.assertEquals(error[0], 0)
        self.assertEquals(len(string_out), 2)
        self.assertEquals(string_out[0], "aaaaa")
        self.assertEquals(string_out[1], "bbbb")
        
        
        del x
        
    def test11(self):
        x = ForTestingInterface()
        string_out, error = x.echo_string(["", "bbbb"])
        self.assertEquals(error[0], 0)
        self.assertEquals(len(string_out), 2)
        self.assertEquals(string_out[0], "")
        self.assertEquals(string_out[1], "bbbb")
        
        
        del x
        
    def test12(self):
        x = ForTestingInterface()
        str1_out, str2_out, error = x.echo_strings("abc", "def")
        self.assertEquals(error, 0)
        self.assertEquals(str1_out[0], "cba")
        self.assertEquals(str2_out[0], "fed")
        del x
        
        
    def test13(self):
        x = ForTestingInterface()
        str1_out, str2_out, error = x.echo_strings(["abc", "def"], ["ghi", "jkl"])
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(str1_out[0], "cba")
        self.assertEquals(str1_out[1], "fed")
        self.assertEquals(str2_out[0], "ihg")
        self.assertEquals(str2_out[1], "lkj")
        del x
        
