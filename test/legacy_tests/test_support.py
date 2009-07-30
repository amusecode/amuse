import unittest
import sys

from amuse.legacy.support.core import *
from amuse.legacy.support.create_c import *
from amuse.legacy.support.create_fortran import *
        

class TestLegacyFunction(unittest.TestCase):
    @legacy_function
    def get_time_step():
        function = RemoteFunction()  
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
    @legacy_function
    def interleave_ints_and_doubles():
        function = RemoteFunction()  
        function.id = 2
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('parameter3', dtype='i', direction=function.INOUT)
        function.addParameter('parameter4', dtype='d', direction=function.INOUT)
        function.addParameter('parameter5', dtype='i', direction=function.IN)
        function.addParameter('parameter6', dtype='d', direction=function.IN)
        return function
        
    def test1(self):
        self.assertTrue(isinstance(self.get_time_step, LegacyCall))
        self.assertTrue(isinstance(type(self).get_time_step, legacy_function))
        
    def test3(self):
        class TestChannel(object):
            def send_message(self, id, doubles_in = [], ints_in = []):
                self.in_doubles = doubles_in
                self.in_ints = ints_in
            
            def recv_message(self, id):
                return ([1,2],[3.0,4.0])
                
        self.channel = TestChannel()
        result = self.get_time_step(1, 2.0, 3.0)
        self.assertFalse(self.channel.in_ints is None)
        self.assertFalse(self.channel.in_doubles is None)
        self.assertEqual(self.channel.in_ints[0] , 1)
        self.assertEqual(self.channel.in_doubles[0] , 2.0)
        self.assertEqual(self.channel.in_doubles[1] , 3.0)
        
    def test4(self):
        class TestChannel(object):
            def send_message(self, id, doubles_in = [], ints_in = []):
                self.in_doubles = doubles_in
                self.in_ints = ints_in
            
            def recv_message(self, id):
                return ([1,2],[3.0,4.0])
                
        self.channel = TestChannel()
        result = self.interleave_ints_and_doubles(1, 2.1, 3, 4.2, 5, 6.3)
        self.assertFalse(self.channel.in_ints is None)
        self.assertFalse(self.channel.in_doubles is None)
        self.assertEqual(self.channel.in_ints[0] , 1)
        self.assertEqual(self.channel.in_ints[1] , 3)
        self.assertEqual(self.channel.in_ints[2] , 5)
        self.assertEqual(self.channel.in_doubles[0] , 2.1)
        self.assertEqual(self.channel.in_doubles[1] , 4.2)
        self.assertEqual(self.channel.in_doubles[2] , 6.3)
        
    def test5(self):
        class TestChannel(object):
            def send_message(self, id, doubles_in = [], ints_in = []):
                self.in_doubles = doubles_in
                self.in_ints = ints_in
            
            def recv_message(self, id):
                return ([1,2],[3.0,4.0])
                
        self.channel = TestChannel()
        result = self.interleave_ints_and_doubles(parameter2 = 2.1, parameter6 = 6.3, parameter1 = 1, parameter4 = 4.2, parameter3 = 3, parameter5 = 5)
        
        self.assertFalse(self.channel.in_ints is None)
        self.assertFalse(self.channel.in_doubles is None)
        
        self.assertEqual(self.channel.in_ints[0] , 1)
        self.assertEqual(self.channel.in_ints[1] , 3)
        self.assertEqual(self.channel.in_ints[2] , 5)
        
        self.assertEqual(self.channel.in_doubles[0] , 2.1)
        self.assertEqual(self.channel.in_doubles[1] , 4.2)
        self.assertEqual(self.channel.in_doubles[2] , 6.3)
        
        
        
class TestMakeACStringOfALegacyFunctionSpecification(unittest.TestCase):
    _class_to_test = MakeACStringOfALegacyFunctionSpecification
    
    def test1(self):
        function = RemoteFunction()      
        function.name = "test_one"     
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        
        x = self._class_to_test()
        
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:test_one(ints_in[0],doubles_in[0],doubles_in[1]);break;')
        
    def test2(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:ints_out[0]=test_one(ints_in[0]);reply.number_of_ints=1;break;')
    
    def test3(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:ints_out[0]=test_one(ints_in[0],&ints_out[1]);reply.number_of_ints=2;break;')
        
    def test4(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.addParameter('doublep', dtype='d', direction=function.OUT)
        function.result_type = 'd'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:doubles_out[0]=test_one(ints_in[0],&ints_out[0],&doubles_out[1]);reply.number_of_ints=1;reply.number_of_doubles=2;break;')
        
    def test6(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.INOUT)
        function.addParameter('doublep', dtype='d', direction=function.INOUT)
        function.result_type = 'd'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:doubles_out[0]=test_one(ints_in[0],&ints_in[1],&doubles_in[0]);ints_out[0]=ints_in[1];doubles_out[1]=doubles_in[0];reply.number_of_ints=1;reply.number_of_doubles=2;break;')
        
    def test5(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.addParameter('doublep', dtype='d', direction=function.OUT)
        function.result_type = 'd'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        self.assertEquals(string,  'case 1:\n  doubles_out[0] = test_one(\n    ints_in[0] ,\n    &ints_out[0] ,\n    &doubles_out[1]\n  );\n  reply.number_of_ints = 1;\n  reply.number_of_doubles = 2;\n  break;')
        
class TestMakeACStringOfAClassWithLegacyFunctions(unittest.TestCase):
    _class_to_test = MakeACStringOfAClassWithLegacyFunctions
    
    @legacy_function
    def get_time_step():
        function = RemoteFunction()  
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
        
    def test1(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeACStringOfAClassWithLegacyFunctions
        string = x.result
        print string
        self.assertTrue('#include <mpi.h>' in string)
        self.assertTrue('#include "parameters.h"' in string)
        self.assertTrue('run_loop' in string)
        
    def test2(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeACStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('main' in string)
        

class TestMakeAFortranStringOfALegacyFunctionSpecification(unittest.TestCase):
    _class_to_test = MakeAFortranStringOfALegacyFunctionSpecification
    
    def test1(self):
        function = RemoteFunction()      
        function.name = "test_one"     
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        
        x = self._class_to_test()
        x.specification = function
        string = x.result
        print string
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(1)CALLtest_one(&integers_in(1),&doubles_in(1),&doubles_in(2)&)')
        
    def test2(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        x = self._class_to_test()
        x.specification = function
        string = x.result
        print string
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(1)CALLtest_one(&integers_in(1),&integers_out(1)&)number_of_integers_out=1')
        
        
        
    def test3(self):
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.addParameter('doublep', dtype='d', direction=function.OUT)
        x = self._class_to_test()
        x.specification = function
        string = x.result
        print string
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(1)\nCALLtest_one(&\nintegers_in(1),&\nintegers_out(1),&\ndoubles_out(1)&\n)\nnumber_of_integers_out=1\nnumber_of_doubles_out=1\n')

    def test4(self):
        
        function = RemoteFunction()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.INOUT)
        function.addParameter('doublep', dtype='d', direction=function.INOUT)
        x = self._class_to_test()
        x.specification = function
        string = x.result
        print string
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(1)\nCALLtest_one(&\nintegers_in(1),&\nintegers_in(2),&\ndoubles_in(1)&\n)\nintegers_out(1)=integers_in(2)\ndoubles_out(1)=doubles_in(1)\nnumber_of_integers_out=1\nnumber_of_doubles_out=1\n')
    
    def xtest5(self):
        from amuse.legacy.sse.muse_stellar_mpi import SSE_new
        function = SSE_new.get_time_step.specification
        
        x = self._class_to_test()
        x.specification = function
        string = x.result
        print string
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(1)\nCALLtest_one(&\nintegers_in(0),&\nintegers_in(1),&\ndoubles_in(0)&\n)\nintegers_out(0)=integers_in(1)\ndoubles_out(0)=doubles_in(0)\nnumber_of_integers_out=1\nnumber_of_doubles_out=1\n')
    
        
      
class TestMakeAFortranStringOfAClassWithLegacyFunctions(unittest.TestCase):
    _class_to_test = MakeAFortranStringOfAClassWithLegacyFunctions
    
    @legacy_function
    def get_time_step():
        function = RemoteFunction()  
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
        
    def test1(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeAFortranStringOfAClassWithLegacyFunctions
        string = x.result
        print string
        self.assertTrue('run_loop' in string)
        
    def test2(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeAFortranStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('program' in string)
          
        
