from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI, TestCase

import numpy
import parser
import sys
import os
import time
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi import python_code
from amuse.rfi.core import *
from amuse.rfi.channel import AsyncRequestsPool

            
class ForTestingInterface(PythonCodeInterface):
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, implementation_factory = ForTestingImplementation, **options)
        
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
        
    @legacy_function
    def sleep():
        function = LegacyFunctionSpecification()  
        function.addParameter('number_of_seconds', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = False
        return function    
    
    @legacy_function
    def sum_doubles():
        function = LegacyFunctionSpecification()
        function.addParameter('double_in1', dtype='float64', direction=function.IN)
        function.addParameter('double_in2', dtype='float64', direction=function.IN)
        function.addParameter('double_out', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
    
    @legacy_function
    def multiply_ints():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in1', dtype='int32', direction=function.IN)
        function.addParameter('int_in2', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
        
    @legacy_function
    def print_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
        
    @legacy_function
    def print_error_string():
        function = LegacyFunctionSpecification()  
        function.addParameter('string_in', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function  
    
    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    
    
    
    
class ForTestingImplementation(object):
    
    def __init__(self):
        self.masses = [0.0] * 100
        self._particle_data = numpy.reshape(numpy.arange(300.0), (-1, 3))
        
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
        
    def print_string(self, string_in):
        print string_in
        return 0
        
    def print_error_string(self, string_in):
        print >> sys.stderr, string_in
        return 0
        
    def echo_strings(self, string_inout1, string_inout2):
        string_inout1.value = string_inout1.value[::-1]
        string_inout2.value = string_inout2.value[::-1]
        return 0
        
    def sleep(self, number_of_seconds):
        import time
        time.sleep(number_of_seconds)
        return 0
    
    def sum_doubles(self, double_in1, double_in2, double_out):
        double_out.value = double_in1 + double_in2
        return 0
    
    def multiply_ints(self, int_in1, int_in2, int_out):
        int_out.value = int_in1 * int_in2
        return 0
    
    def get_position(self, index_of_the_particle, x, y, z, length):
        try:
            x.value = self._particle_data[index_of_the_particle, 0]
            y.value = self._particle_data[index_of_the_particle, 1]
            z.value = self._particle_data[index_of_the_particle, 2]
            return 0
        except:        
            return -1
    
    def set_position(self, index_of_the_particle, x, y, z, length):
        try:
            self._particle_data[index_of_the_particle, 0] = x
            self._particle_data[index_of_the_particle, 1] = y
            self._particle_data[index_of_the_particle, 2] = z
            return 0
        except:        
            return -1
    

class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(**options), **options)
    
    def define_methods(self, object):
        object.add_method("sleep", (units.s,), (object.ERROR_CODE,))


            
class TestInterface(TestWithMPI):
    def setUp(self):
        self.check_not_in_mpiexec()
        super(TestInterface, self).setUp()
        
    def check_not_in_mpiexec(self):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 
        
        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            self.skip('cannot run the socket tests under hydra process manager')
            
            
    def test2(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.SocketMessage(0, 10, 1)
        input_message.ints = [1]
        
        output_message = python_code.SocketMessage(0, 10, 1)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 1)
        self.assertEquals(len(output_message.doubles), 1)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(output_message.doubles[0], 0.0)
        
    def test3(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.SocketMessage(0, 11, 1)
        input_message.ints = [1]
        input_message.doubles = [12.0]
        
        output_message = python_code.SocketMessage(0, 10, 1)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 1)
        self.assertEquals(len(output_message.doubles), 0)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(implementation.masses[1], 12.0)
        
    
    def test4(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.SocketMessage(0, 11, 4)
        input_message.ints = [1,2,3,4]
        input_message.doubles = [12.0,13.0,14.0,15.0]
        
        output_message = python_code.SocketMessage(0, 10, 4)
        
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
        x = ForTestingInterface(redirection='none', channel_type = 'sockets')
        
        error = x.set_mass(1, 10.0)
        self.assertEquals(error, 0)
        
        answer, error = x.get_mass(1)
        self.assertEquals(error, 0)
        self.assertEquals(answer, 10.0)
        x.stop()
        
        
        
    def test6(self):
        x = ForTestingInterface(channel_type = 'sockets')
        
        errors = x.set_mass([1,2], [10.0,11.0])
        self.assertEquals(errors[0], 0)
        self.assertEquals(errors[1], 0)
        
        answer, errors = x.get_mass([1,2])
        self.assertEquals(errors[0], 0)
        self.assertEquals(answer[0], 10.0)
        self.assertEquals(answer[1], 11.0)
        x.stop()
        
        x.stop()
        
    def test7(self):
        x = ForTestingInterface(channel_type = 'sockets')
        
        int_out, error = x.echo_int(20)
        self.assertEquals(error, 0)
        self.assertEquals(int_out, 20)
        x.stop()
        
        x.stop()
        
    
    def test8(self):
        implementation = ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, ForTestingInterface)
        
        input_message = python_code.SocketMessage(0, 12, 1)
        input_message.ints = [20]
        
        output_message = python_code.SocketMessage(0, 10, 1)
        
        x.handle_message(input_message, output_message)
        
        self.assertEquals(len(output_message.ints), 2)
        self.assertEquals(output_message.ints[0], 0)
        self.assertEquals(output_message.ints[1], 20)
        
    
    def test9(self):
        x = ForTestingInterface(channel_type = 'sockets')
        string_out, error = x.echo_string("1234567")
        self.assertEquals(error, 0)
        self.assertEquals(string_out, "1234567")
        x.stop()
        
    def test10(self):
        x = ForTestingInterface(channel_type = 'sockets')
        string_out, error = x.echo_string(["aaaaa", "bbbb"])
        self.assertEquals(error[0], 0)
        self.assertEquals(len(string_out), 2)
        self.assertEquals(string_out[0], "aaaaa")
        self.assertEquals(string_out[1], "bbbb")
        x.stop()
        
    def test11(self):
        x = ForTestingInterface(channel_type = 'sockets')
        string_out, error = x.echo_string(["", "bbbb"])
        self.assertEquals(error[0], 0)
        self.assertEquals(len(string_out), 2)
        self.assertEquals(string_out[0], "")
        self.assertEquals(string_out[1], "bbbb")
        x.stop()
        
    def test12(self):
        x = ForTestingInterface(channel_type = 'sockets')
        str1_out, str2_out, error = x.echo_strings("abc", "def")
        self.assertEquals(error, 0)
        self.assertEquals(str1_out, "cba")
        self.assertEquals(str2_out, "fed")
        x.stop()
        
        
    def test13(self):
        x = ForTestingInterface(channel_type = 'sockets')
        str1_out, str2_out, error = x.echo_strings(["abc", "def"], ["ghi", "jkl"])
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(str1_out[0], "cba")
        self.assertEquals(str1_out[1], "fed")
        self.assertEquals(str2_out[0], "ihg")
        self.assertEquals(str2_out[1], "lkj")
        x.stop()
        
    def test14(self):
        x = ForTestingInterface(channel_type = 'sockets')
        result = x.sleep(0.01)
        self.assertEquals(result, 0)
        request = x.sleep.async(0.01)
        request.wait()
        result = request.result()
        self.assertEquals(result, 0)
        x.stop()
        
    def test15(self):
        x = ForTestingInterface(channel_type = 'sockets')
        y = ForTestingInterface(channel_type = 'sockets')
        request1 = x.sleep.async(0.5)
        self.assertFalse(request1.is_result_available())
        request2 = y.sleep.async(1.5)
        self.assertFalse(request2.is_result_available())
        request2.wait()
        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())
        
        self.assertEquals(request1.result(), 0)
        self.assertEquals(request2.result(), 0)
        
        y.stop()
        x.stop()
        
    
    
    def test16(self):
        x = ForTestingInterface(channel_type = 'sockets')
        request1 = x.sleep.async(0.4)
        self.assertRaises(Exception, lambda : x.sleep(0.01))
        request1.wait()
        self.assertRaises(Exception, lambda : x.sleep(0.01))
        request1.result()
        x.sleep(0.01)
        self.assertTrue(request1.is_result_available())
        x.stop()
    
    def test17(self):
        x = ForTesting(channel_type = 'sockets')
        self.assertTrue(x.sleep.is_async_supported)
        request= x.sleep.async(0.2 | units.s)
        request.wait()
        result = request.result()
        self.assertEquals(result, [])
        x.stop()
    
    
    
    def test21(self):
        print "Testing must_handle_array for Python codes"
        instance = ForTestingInterface(channel_type = 'sockets')
        
        x,y,z,err = instance.get_position(range(100))
        self.assertEquals(err, 0)
        self.assertEquals(x, numpy.arange(0.0, 300.0, 3.0))
        self.assertEquals(y, numpy.arange(1.0, 300.0, 3.0))
        self.assertEquals(z, numpy.arange(2.0, 300.0, 3.0))
        x,y,z,err = instance.get_position(range(101))
        self.assertEquals(err, -1)
        self.assertEquals(instance.get_position(1).values(), [3.0, 4.0, 5.0, 0])
        
        err = instance.set_position(range(100), numpy.arange(100.0), numpy.arange(100.0, 200.0), numpy.arange(200.0, 300.0))
        self.assertEquals(err, 0)
        err = instance.set_position(range(101), numpy.arange(101.0), numpy.arange(101.0, 202.0), numpy.arange(202.0, 303.0))
        self.assertEquals(err, -1)
        err = instance.set_position(0, -1.0, -2.0, -3.0)
        
        x,y,z,err = instance.get_position(range(100))
        self.assertEquals(err, 0)
        self.assertEquals(x, numpy.concatenate(([-1.0], numpy.arange(1.0, 100.0))))
        self.assertEquals(y, numpy.concatenate(([-2.0], numpy.arange(101.0, 200.0))))
        self.assertEquals(z, numpy.concatenate(([-3.0], numpy.arange(201.0, 300.0))))

        instance.stop()
    
    
    def test22(self):
        
        pool = AsyncRequestsPool()
        
        x = ForTestingInterface(channel_type = 'sockets')
        y = ForTestingInterface(channel_type = 'sockets')
        request1 = x.sleep.async(0.5)
        request2 = y.sleep.async(1.5)
        finished_requests = []
        
        def handle_result(request, index):
            self.assertTrue(request.is_result_available())
            finished_requests.append(index)
            
        pool.add_request(request1, handle_result, [1])
        pool.add_request(request2, handle_result, [2])
        
        pool.wait()
        self.assertEquals(len(finished_requests), 1)
        self.assertEquals(len(pool), 1)
        
        pool.wait()
        self.assertEquals(len(finished_requests), 2)
        self.assertEquals(len(pool), 0)
        
        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())
        
        self.assertEquals(request1.result(), 0)
        self.assertEquals(request2.result(), 0)
        
        y.stop()
        x.stop()

    def test23(self):
        
        pool = AsyncRequestsPool()
        
        x = ForTestingInterface(channel_type = 'sockets')
        y = ForTestingInterface(channel_type = 'sockets')
        request1 = x.sleep.async(0.2)
        request2 = y.sleep.async(0.2)
        finished_requests = []
        
        def handle_result(request, index):
            self.assertTrue(request.is_result_available())
            finished_requests.append(index)
            
        pool.add_request(request1, handle_result, [1])
        pool.add_request(request2, handle_result, [2])
        
        time.sleep(1.0)

        pool.wait()
        pool.wait()

        self.assertEquals(len(finished_requests), 2)
        self.assertEquals(len(pool), 0)
        
        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())
        
        self.assertEquals(request1.result(), 0)
        self.assertEquals(request2.result(), 0)
        
        pool.wait()
        self.assertEquals(len(pool), 0)

        y.stop()
        x.stop()



    def test18(self):
        print "Testing the splitting of very long MPI messages into blocks"
        x = ForTesting(max_message_length=10, channel_type = 'sockets')
        N = 100
        doubles = x.echo_double([1.0*i for i in range(N)])
        self.assertTrue(list(doubles) == [1.0*i for i in range(N)])
        sums = x.sum_doubles([1.0*i for i in range(N)],[1.0*i for i in range(N)])
        self.assertTrue(list(sums) == [2.0*i for i in range(N)])
        products = x.multiply_ints(range(N),range(N))
        self.assertTrue(list(products) == [i*i for i in range(N)])
        N = 101
        doubles = x.echo_double([1.0*i for i in range(N)])
        self.assertTrue(list(doubles) == [1.0*i for i in range(N)])
        sums = x.sum_doubles([1.0*i for i in range(N)],[1.0*i for i in range(N)])
        self.assertTrue(list(sums) == [2.0*i for i in range(N)])
        products = x.multiply_ints(range(N),range(N))
        self.assertTrue(list(products) == [i*i for i in range(N)])
        x.stop()
        
        


