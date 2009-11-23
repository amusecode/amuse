import unittest
import sys

from amuse.legacy.support.core import *
from amuse.legacy.support.create_c import *
from amuse.legacy.support.create_fortran import *
        

class TestLegacyFunction(unittest.TestCase):
    @legacy_function
    def get_time_step():
        function = LegacyFunctionSpecification()  
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
    @legacy_function
    def interleave_ints_and_doubles():
        function = LegacyFunctionSpecification()  
        function.id = 2
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('parameter3', dtype='i', direction=function.INOUT)
        function.addParameter('parameter4', dtype='d', direction=function.INOUT)
        function.addParameter('parameter5', dtype='i', direction=function.IN)
        function.addParameter('parameter6', dtype='d', direction=function.IN)
        return function
    
    @legacy_function
    def send_string():
        function = LegacyFunctionSpecification()
        function.id = 3
        function.addParameter('parameter1', dtype='string', direction=function.IN)
        return function
        
    def test1(self):
        self.assertTrue(isinstance(self.get_time_step, LegacyCall))
        self.assertTrue(isinstance(type(self).get_time_step, legacy_function))
        
    def test3(self):
        class TestChannel(object):
            def send_message(self, id, doubles_in = [], ints_in = []):
                self.in_doubles = doubles_in
                self.in_ints = ints_in
            
            def recv_message(self, id, handle_as_array):
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
            
            def recv_message(self, id, handle_as_array):
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
            
            def recv_message(self, id, handle_as_array):
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
        
    
    def test6(self):
        class TestChannel(object):
            def send_message(self, id, doubles_in = [], ints_in = [], chars_in=[]):
                self.chars_in = chars_in
            
            def recv_message(self, id, handle_as_array):
                return [[], []]
                
        self.channel = TestChannel()
        result = self.send_string('bla')
        
        self.assertFalse(self.channel.chars_in is None)
        self.assertEqual(self.channel.chars_in[0], 'bla')
        
        
        
        
        
class TestMakeACStringOfALegacyFunctionSpecification(unittest.TestCase):
    _class_to_test = MakeACStringOfALegacyFunctionSpecification
    
    def test1(self):
        function = LegacyFunctionSpecification()      
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
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:ints_out[0]=test_one(ints_in[0]);reply_header.number_of_ints=1;break;')
    
    def test3(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:ints_out[0]=test_one(ints_in[0],&ints_out[1]);reply_header.number_of_ints=2;break;')
    
   
    def test4(self):
        function = LegacyFunctionSpecification()      
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
        self.assertEquals(string_no_spaces, 'case1:doubles_out[0]=test_one(ints_in[0],&ints_out[0],&doubles_out[1]);reply_header.number_of_ints=1;reply_header.number_of_doubles=2;break;')
    
  
    def test5(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.addParameter('doublep', dtype='d', direction=function.OUT)
        function.result_type = 'd'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        self.assertEquals(string,  'case 1:\n  doubles_out[0] = test_one(\n    ints_in[0] ,\n    &ints_out[0] ,\n    &doubles_out[1]\n  );\n  reply_header.number_of_ints = 1;\n  reply_header.number_of_doubles = 2;\n  break;')


    def test6(self):
        function = LegacyFunctionSpecification()      
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
        self.assertEquals(string_no_spaces, 'case1:doubles_out[0]=test_one(ints_in[0],&ints_in[1],&doubles_in[0]);ints_out[0]=ints_in[1];doubles_out[1]=doubles_in[0];reply_header.number_of_ints=1;reply_header.number_of_doubles=2;break;')
           
    def test7(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='string', direction=function.IN)
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case1:test_one(characters+(0-1<0?0:strings_in[0-1]+1));break;')
        

    def test8(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.INOUT)
        function.addParameter('doublep', dtype='d', direction=function.INOUT)
        function.result_type = 'd'
        function.can_handle_array = True
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        expected = """
        case 1:
          for (int i = 0 ; i < request_header.len; i++){
            doubles_out[i] = test_one(
              ints_in[i] ,
              &ints_in[( 1 * request_header.len) + i] ,
              &doubles_in[i]
            );
            ints_out[i] = ints_in[( 1 * request_header.len) + i];
            doubles_out[( 1 * request_header.len) + i] = doubles_in[i];
          }
          reply_header.number_of_ints = 1;
          reply_header.number_of_doubles = 2;
          break;
        """
        print string
        expected_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , expected))
        self.assertEquals(string_no_spaces, expected_no_spaces);
        

    def test9(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.can_handle_array = True
        function.addParameter('parameter1', dtype='string', direction=function.IN)
        function.addParameter('parameter2', dtype='string', direction=function.IN)
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        expected = """
        case 1:
          for (int i = 0 ; i < request_header.len; i++){
            test_one(
              characters + ( i- 1 < 0 ? 0 :strings_in[i - 1] + 1) ,
              characters + ( ( 1 * request_header.len) + i- 1 < 0 ? 0 :strings_in[( 1 * request_header.len) + i - 1] + 1)
            );
          }
          break;
        """
        print string
        expected_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , expected))
        self.assertEquals(string_no_spaces, expected_no_spaces);

        

        
class TestMakeACHeaderDefinitionStringOfALegacyFunctionSpecification(unittest.TestCase):
    _class_to_test = MakeACHeaderDefinitionStringOfALegacyFunctionSpecification
    
    def test1(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"     
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        
        x = self._class_to_test()
        
        x.specification = function
        string = x.result
        self.assertEquals(string, 'void test_one(int parameter1, double parameter2, double name);')
        
    def test2(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.result_type = 'i'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        self.assertEquals(string, 'int test_one(int parameter1);')
    
    def test3(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_with_out"   
        function.id = 1
        function.addParameter('parameter1', dtype='int32', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.OUT)
        function.addParameter('doublep', dtype='float64', direction=function.OUT)
        function.result_type = 'd'
        x = self._class_to_test()
        x.specification = function
        string = x.result
        self.assertEquals(string,  'double test_with_out(int parameter1, int * parameter2, double * doublep);')


        
class TestMakeACStringOfALegacyGlobalSpecification(unittest.TestCase):
    _class_to_test = MakeACStringOfALegacyGlobalSpecification
    
    def test1(self):
        x = self._class_to_test()
        x.legacy_global = legacy_global('test',id=10,dtype='d')
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case10:if(request_header.number_of_doubles==1){test=doubles_in[0];}else{reply_header.number_of_doubles=1;doubles_out[0]=test;}break;')
    
    def test2(self):
        x = self._class_to_test()
        x.legacy_global = legacy_global('test',id=10,dtype='i')
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        self.assertEquals(string_no_spaces, 'case10:if(request_header.number_of_ints==1){test=ints_in[0];}else{reply_header.number_of_ints=1;ints_out[0]=test;}break;')
        


class TestMakeACStringOfAClassWithLegacyFunctions(unittest.TestCase):
    _class_to_test = MakeACStringOfAClassWithLegacyFunctions
    
    include_headers = ['muse_dynamics.h', 'parameters.h', 'local.h']
    
    extra_content = """
    int _add_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz)) {
        dynamics_state state;
        state.id = id;
        state.mass = mass;
        state.radius = radius;
        state.x = x;
        state.y = y
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        return add_particle(state);
    }
    
    void _get_state(int id, int * id_out,  double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz) {
        dynamics_state state = get_state(id)
        *id_out = state.id
        *mass = state.mass;
        *radius = state.radius;
        *x = state.x;
        *y = state.y;
        *z = state.z;
        *vx = state.vx;
        *vy = state.vy;
        *vz = state.vz;
    }
    """
    
    
    @legacy_function
    def get_time_step():
        function = LegacyFunctionSpecification()  
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
        
    def test1(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeACStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('#include <mpi.h>' in string)
        self.assertTrue('#include "parameters.h"' in string)
        self.assertTrue('run_loop' in string)
        
    def test2(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeACStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('main' in string)
    
    def test3(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeACStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('_get_state(' in string)
        
    def test4(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeACStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('int number_of_ints;' in string)
        self.assertTrue('number_of_ints(0)' in string)
        self.assertTrue('int header[6];' in string)

class TestMakeAFortranStringOfALegacyFunctionSpecification(unittest.TestCase):
    _class_to_test = MakeAFortranStringOfALegacyFunctionSpecification
    
    def test1(self):
        function = LegacyFunctionSpecification()      
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
        function = LegacyFunctionSpecification()      
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
        function = LegacyFunctionSpecification()      
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
        
        function = LegacyFunctionSpecification()      
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
        from amuse.legacy.sse.muse_stellar_mpi import SSE
        function = SSE.evolve.specification
        
        x = self._class_to_test()
        x.specification = function
        string = x.result
        print string
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(1)\nCALLtest_one(&\nintegers_in(0),&\nintegers_in(1),&\ndoubles_in(0)&\n)\nintegers_out(0)=integers_in(1)\ndoubles_out(0)=doubles_in(0)\nnumber_of_integers_out=1\nnumber_of_doubles_out=1\n')
    
    def test6(self):
        function = LegacyFunctionSpecification()      
        function.name = "test_one"   
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='i', direction=function.INOUT)
        function.addParameter('doublep', dtype='d', direction=function.INOUT)
        function.result_type = 'd'
        function.can_handle_array = True
        x = self._class_to_test()
        x.specification = function
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , string))
        expected = """
        CASE(1)
          do i = 1, len_in, 1
            doubles_out(i) = test_one( &
              integers_in(i) ,&
              integers_in(( 1 * len_in) + i) ,&
              doubles_in(i) &
            )
            integers_out(i) = integers_in(( 1 * len_in) + i)
            doubles_out(( 1 * len_in) + i) = doubles_in(i)
          end do
          number_of_integers_out = 1
          number_of_doubles_out = 2
        """
        print string
        expected_no_spaces = ''.join(filter(lambda x : x not in ' \t\n\r' , expected))
        self.assertEquals(string_no_spaces, expected_no_spaces);
      
class TestMakeAFortranStringOfAClassWithLegacyFunctions(unittest.TestCase):
    _class_to_test = MakeAFortranStringOfAClassWithLegacyFunctions
    
    @legacy_function
    def get_time_step():
        function = LegacyFunctionSpecification()  
        function.id = 1
        function.addParameter('parameter1', dtype='i', direction=function.IN)
        function.addParameter('parameter2', dtype='d', direction=function.IN)
        function.addParameter('name', dtype='d', direction=function.IN)
        return function
        
    
    @legacy_function
    def get_number():
        function = LegacyFunctionSpecification()  
        function.id = 2
        function.result_type = 'i'
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
        
    def test3(self):
        x = self._class_to_test()
        x.class_with_legacy_functions = TestMakeAFortranStringOfAClassWithLegacyFunctions
        string = x.result
        self.assertTrue('integer :: get_number' in string)
        
                        
class TestMakeAFortranStringOfALegacyGlobalSpecification(unittest.TestCase):
    _class_to_test = MakeAFortranStringOfALegacyGlobalSpecification
    
    def test1(self):
        x = self._class_to_test()
        x.legacy_global = legacy_global('test',id=10,dtype='d')
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(10)\nif(number_of_doubles_in==1)then\ntest=doubles_in[1]\nelse\nnumber_of_doubles_out=1\ndoubles_out[1]=test\n\n')
    
    def test2(self):
        x = self._class_to_test()
        x.legacy_global = legacy_global('test',id=10,dtype='i')
        string = x.result
        string_no_spaces = ''.join(filter(lambda x : x not in ' \t\r' , string))
        self.assertEquals(string_no_spaces, 'CASE(10)\nif(number_of_integers_in==1)then\ntest=integers_in[1]\nelse\nnumber_of_integers_out=1\nintegers_out[1]=test\n\n')
        


          
        
