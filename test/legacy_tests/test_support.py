import unittest
import sys

from amuse.legacy.support.core import *
        

class TestLegacyFunction(unittest.TestCase):
    @legacy_function
    def get_time_step():
        function = RemoteFunction()      
        function.name = "test_one"
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
        
    def test1(self):
        self.assertTrue(isinstance(self.get_time_step, legacy_call))
        self.assertTrue(isinstance(type(self).get_time_step, legacy_function))
        
    def test2(self):
        f = type(self).get_time_step
        string = f.to_c_string()
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'test_one(ints_in[0],doubles_in[0],doubles_in[1]);')
        
class TestMakeACStringOfALegacyFunctionSpecification(unittest.TestCase):
    _class_to_test = MakeACStringOfALegacyFunctionSpecification
    
    def test1(self):
        function = RemoteFunction()      
        function.name = "test_one"
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'test_one(ints_in[0],doubles_in[0],doubles_in[1]);')
        
    def test2(self):
        function = RemoteFunction()      
        function.name = "test_one"
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'ints_out[0]=test_one(ints_in[0]);')
    
    def test3(self):
        function = RemoteFunction()      
        function.name = "test_one"
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'ints_out[0]=test_one(ints_in[0],&ints_out[1]);')
        
        
        
        
        