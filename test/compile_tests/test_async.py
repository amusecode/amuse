from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions
from amuse.support import options

import os
import time
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi.tools import create_c
from amuse.rfi import async_request
from amuse.rfi.core import *

from . import test_c_implementation

from amuse.test import compile_tools

from amuse.community.distributed.interface import DistributedAmuse, Pilot

codestring=test_c_implementation.codestring+"""
#include <unistd.h>

float _x[10] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

int do_sleep(int in) {
    sleep(in);
    return 0;
}

int return_error(int * out) {
    *out=123;
    return -1;
}

int get_x(int in, float *x){
    *x=_x[in];
    return 0;
    }

int set_x(int in, float x){
    _x[in]=x;
    return 0;
    }

int dummy(){
    return 0;
    }

"""

class ForTestingInterface(test_c_implementation.ForTestingInterface):
    @legacy_function
    def do_sleep():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function 

    @legacy_function
    def return_error():
        function = LegacyFunctionSpecification()
        function.addParameter('out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function 

    @legacy_function
    def echo_2_int():
        function = LegacyFunctionSpecification()
        function.addParameter('int_in1', dtype='int32', direction=function.IN, unit=units.m)
        function.addParameter('int_in2', dtype='int32', direction=function.IN, default = 1, unit=units.kg)
        function.addParameter('int_out1', dtype='int32', direction=function.OUT, unit=units.m)
        function.addParameter('int_out2', dtype='int32', direction=function.OUT, unit=units.kg)
        function.addParameter('len', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def get_x():
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float32', direction=function.OUT, unit=units.m)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function 

    @legacy_function
    def set_x():
        function = LegacyFunctionSpecification()
        function.addParameter('index', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float32', direction=function.IN, unit=units.m)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function 

    @legacy_function
    def dummy():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function 

class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)

    def get_grid_range(self):
        return (0,9)

    def define_grids(self, handler):
        handler.define_grid('grid')
        handler.set_grid_range('grid', 'get_grid_range')
        handler.add_getter('grid', 'get_x', names=["x"])
        handler.add_setter('grid', 'set_x', names=["x"])

class ForTestingWithState(ForTesting):
    def define_state(self, handler):
        handler.set_initial_state("1")
        handler.add_transition("1", "2", "dummy")
        handler.add_method("2", "get_x")
        handler.add_method("2", "set_x")

class TestASync(TestWithMPI):

    @classmethod
    def setup_class(cls):
        print("building...")
        cls.check_can_compile_modules()
        try:
            cls.exefile = compile_tools.build_worker(codestring, cls.get_path_to_results(), ForTestingInterface)
        except Exception as ex:
            print(ex)
            raise
        print("done")
        
    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEqual(int_out, 10)
        self.assertEqual(error, 0)

    def test2(self):
        instance = ForTestingInterface(self.exefile)
        request = instance.echo_int.asynchronous(10)
        self.assertEqual(request, instance.async_request)
        request.wait()
        int_out,error=request.result()
        self.assertEqual(int_out, 10)
        self.assertEqual(error, 0)
        instance.stop()

    def test3(self):
        instance = ForTestingInterface(self.exefile)
        request1 = instance.do_sleep.asynchronous(1)
        request2 = instance.echo_int.asynchronous(10)
        self.assertEqual(request2, instance.async_request)
        request2.wait()
        int_out,error=request2.result()
        self.assertEqual(int_out, 10)
        self.assertEqual(error, 0)
        instance.stop()

    def test4(self):
        instance = ForTesting(self.exefile)
        request1 = instance.do_sleep(1, return_request=True)
        request2 = instance.echo_int(10, return_request=True)
        self.assertEqual(request2, instance.async_request)
        instance.async_request.wait()
        int_out=request2.result()
        self.assertEqual(int_out, 10)
        instance.stop()

    def test5(self):
        instance = ForTesting(self.exefile)
        instance.do_sleep(1, return_request=True)
        requests=[]
        for x in range(10):
            requests.append(instance.echo_int(x, return_request=True))
        instance.async_request.wait()
        for i,x in enumerate(requests):
            self.assertEqual(x.result(), i)
        instance.stop()

    def test6(self):
        instance = ForTesting(self.exefile)
        requests=[]
        for x in range(10):
            requests.append(instance.echo_int(x, return_request=True))
        instance.async_request.wait()
        for i,x in enumerate(requests):
            self.assertEqual(x.result(), i)
        instance.stop()

    def test7(self):
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        t1=time.time()

        requests=[]
        for x in range(10):
            requests.append([instance1.echo_int(x, return_request=True),x])
        for x in range(10):
            requests.append([instance2.echo_int(x, return_request=True),x])

        instance1.do_sleep(1, return_request=True)
        instance2.do_sleep(1, return_request=True)

        pool=instance1.async_request.join(instance2.async_request)
        pool.waitall()
        t2=time.time()

        for x in requests:
            self.assertEqual(x[0].result(), x[1])
        instance1.stop()
        instance2.stop()
        self.assertTrue(t2-t1 < 2.)

    def test8(self):
        from threading import Thread
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        t1=time.time()

        requests=[]
        for x in range(10):
            requests.append([instance1.echo_int(x, return_request=True),x])
        for x in range(10):
            requests.append([instance2.echo_int(x, return_request=True),x])

        instance1.do_sleep(1, return_request=True)
        instance2.do_sleep(1, return_request=True)

        pool=instance1.async_request.join(instance2.async_request)
        
        thread=Thread(target=pool.waitall)
        thread.start()
        time.sleep(1)
        thread.join()
        
        self.assertTrue(pool)
        
        t2=time.time()

        for x in requests:
            self.assertEqual(x[0].result(), x[1])
        instance1.stop()
        instance2.stop()
        self.assertTrue(t2-t1 < 2.)

    def test9(self):
        instance = ForTesting(self.exefile)
        for x in range(10):
            instance.echo_int(x, return_request=True)
        results=instance.async_request.results
        self.assertEqual(results, list(range(10)))
        instance.stop()

    def test10(self):
        instance = ForTesting(self.exefile)
        r1=instance.do_sleep(1, return_request=True)
        r2=instance.return_error( return_request=True)
        r3=instance.echo_int(1, return_request=True)
        r4=instance.echo_int(2, return_request=True)
        self.assertRaises(Exception, instance.async_request.waitall,
            expected_message="Error in dependent call: Error when calling 'return_error' of a 'ForTesting', errorcode is -1")
        self.assertTrue( r1.is_result_available() )
        self.assertFalse( r2.is_result_available() )
        self.assertTrue(  r2.is_finished )
        self.assertTrue(  r3.is_finished )
        self.assertFalse(  bool(r3) )
        self.assertTrue(  r4.is_finished )
        self.assertTrue(  r4.waits_for() is None )
        self.assertFalse( r3.is_result_available() )
        
        instance.stop()

    def test11(self):
        """ cross dependency """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_int(10, return_request=True)
        
        def fac():
          return instance2.echo_int(20, return_request=True)

        #~ request2=instance2.echo_int(20, async_dependency=request1)
        request2=async_request.DependentASyncRequest(request1, fac)

        request2.wait()
        
        self.assertEqual(request2.result(),20)
        
        instance1.stop()
        instance2.stop()

    def test11b(self):
        """ cross dependency """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_int(10, return_request=True)
        request2=instance2.echo_int(20, async_dependency=request1, return_request=True)

        request2.wait()
        self.assertTrue(request1.is_result_available())
        
        self.assertEqual(request2.result(),20)
        
        instance1.stop()
        instance2.stop()

    def test12(self):
        """ cross dependency with input-output dependency """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_int(10, return_request=True)
        
        results=dict()
        
        def safe_result(arg, index):
            result=arg()
            results[index]=result
            return result
        
        request1.add_result_handler(safe_result,(1,))
        
        def fac():
          return instance2.echo_int(results[1], return_request=True)

        #~ request2=instance2.echo_int(??, async_factory=fac)
        request2=async_request.DependentASyncRequest(request1, fac)

        request2.wait()
        
        self.assertEqual( request2.result(), 10)
        
        instance1.stop()
        instance2.stop()

    def test12b(self):
        """ cross dependency with input-output dependency """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_int(10, return_request=True)

        request2=instance2.echo_int(request1, return_request=True)
        
        request2.wait()

        self.assertEqual( request2.result(), 10)
        
        instance1.stop()
        instance2.stop()

    def test12c(self):
        """ cross dependency with input-output dependency """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_2_int(1 | units.m , 2 | units.kg, return_request=True)

        request2=instance2.echo_2_int(request1[0], request1[1], return_request=True)
        
        print("do...wait...")
        request2.wait()
        print("done", request2.result())

        self.assertEqual( request2.result()[0], 1 | units.m)
        self.assertEqual( request2.result()[1], 2 | units.kg)
        
        instance1.stop()
        instance2.stop()

    def test12c(self):
        """ cross dependency with input-output dependency """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance3 = ForTesting(self.exefile)
        instance4 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_2_int(1 | units.m , 2 | units.kg, return_request=True)
        request1b=instance1.do_sleep(1, return_request=True)
        request2=instance2.echo_2_int(3 | units.m , 4 | units.kg, return_request=True)
        request3=instance3.echo_2_int(request2[0] , 5 | units.kg, return_request=True)

        instance4.do_sleep(1, return_request=True)
        request4=instance4.echo_2_int(request2[0], request3[1], return_request=True, async_dependency=request1b)
        
        request3.wait()

        self.assertEqual( request4.result()[0], 3 | units.m)
        self.assertEqual( request4.result()[1], 5 | units.kg)
        
        instance1.stop()
        instance2.stop()
        instance3.stop()
        instance4.stop()


    def test13(self):
        instance = ForTesting(self.exefile)

        r=instance.echo_int(1, return_request=True)
        time.sleep(0.1)
        self.assertTrue(r.is_result_available())
        r.result()

        r=instance.return_error(return_request=True)
        time.sleep(0.1)
        self.assertTrue(r.is_result_available())
        self.assertTrue(r.is_result_available())
        self.assertRaises(Exception, r.result, expected_message="Error when calling 'return_error' of a 'ForTesting', errorcode is -1")        
        self.assertFalse(r.is_result_available())
        self.assertTrue(r.is_finished)
        instance.stop()

    def test14(self):
        instance = ForTesting(self.exefile, channel_type="sockets")

        r=instance.echo_int(1, return_request=True)
        time.sleep(0.1)
        self.assertTrue(r.is_result_available())
        r.result()

        r=instance.return_error(return_request=True)
        time.sleep(0.1)
        self.assertTrue(r.is_result_available())
        self.assertTrue(r.is_result_available())
        self.assertRaises(Exception, r.result, expected_message="Error when calling 'return_error' of a 'ForTesting', errorcode is -1")        
        self.assertFalse(r.is_result_available())
        self.assertTrue(r.is_finished)
        instance.stop()

    def test15(self):
        instance = ForTesting(self.exefile)

        instance.do_sleep(1, return_request=True)
        instance.return_error( return_request=True)
        instance.echo_int(1, return_request=True)
        instance.echo_int(1, return_request=True)
        #~ self.assertRinstance.echo_int(1)
        

        instance.stop()

    def test16(self):
        instance = ForTesting(self.exefile)

        instance.do_sleep(1, return_request=True)
        result=instance.echo_2_int([11,12,13] | units.m,[3,2,1]| units.kg, return_request=True)

        r1=result[0]
        r2=result[1]
        
        self.assertEqual(r1.result(),[11,12,13] | units.m)
        self.assertEqual(r2.result(),[3,2,1] | units.kg)
        
        instance.stop()

    def test17(self):
        instance = ForTestingInterface(self.exefile)

        instance.do_sleep.asynchronous(1)
        request=instance.echo_2_int.asynchronous([11,12,13],[3,2,1])

        r1=request["int_out1"]
        r2=request["int_out2"]
        
        self.assertEqual(r1.result(),[11,12,13] )
        self.assertEqual(r2.result(),[3,2,1] )
        

        instance.stop()

    def test18(self):
        """ test pool as depedency 1 """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance3 = ForTesting(self.exefile)
        
        request0=instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_int(10, return_request=True)
        request2=instance2.echo_int(10, return_request=True)
        request=async_request.AsyncRequestsPool(request1,request2)
                        
        request3=instance3.echo_int(11, async_dependency=request, return_request=True)
        request3.wait()
        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())
        self.assertEqual( request3.result(), 11)
        
        instance1.stop()
        instance2.stop()
        instance3.stop()

    def test18b(self):
        """ test pool as depedency 2 """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance3 = ForTesting(self.exefile)
        
        request0=instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_int(10, return_request=True)
        request2=instance1.echo_int(10, return_request=True)
        request=async_request.AsyncRequestsPool(request1,request2)
                        
        request3=instance3.echo_int(11, async_dependency=request, return_request=True)
        request3.wait()
        self.assertTrue(request1.is_result_available())
        self.assertTrue(request2.is_result_available())
        self.assertEqual( request3.result(), 11)
        
        instance1.stop()
        instance2.stop()
        instance3.stop()

    def test19(self):
        """ test sum request """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        
        r1=instance1.echo_int(1, return_request=True)
        r2=instance2.echo_int(2, return_request=True)
        s=r1+r2
        r1=instance1.echo_int(2, return_request=True)
        r2=instance2.echo_int(3, return_request=True)
        m=r1*r2
        r1=instance1.echo_int(12, return_request=True)
        r2=instance2.echo_int(3, return_request=True)
        d=r1/r2
        
        self.assertEqual( s.result(), 3)
        self.assertEqual( m.result(), 6)
        self.assertEqual( d.result(), 4)
        
        instance1.stop()
        instance2.stop()

    def test19b(self):
        """ test sum request """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        
        r1=instance1.echo_int(1, return_request=True)
        s=r1+2
        r1=instance1.echo_int(2, return_request=True)
        m=r1*3
        r2=instance2.echo_int(3, return_request=True)
        d=12/r2
        
        self.assertEqual( s.result(), 3)
        self.assertEqual( m.result(), 6)
        self.assertEqual( d.result(), 4)
        
        instance1.stop()
        instance2.stop()

    def test20(self):
        """ some more tests of request expressions """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        instance3 = ForTesting(self.exefile)
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_2_int(1 | units.m , 2 | units.kg, return_request=True)
        request2=instance2.echo_2_int(4 | units.m , 6 | units.kg, return_request=True)
        request3a=instance3.echo_int(request2[0] / request1[0]-4, return_request=True)
        request3b=instance3.echo_int(request2[1] / request1[1]-3, return_request=True)

        request3a.wait()
        request3b.wait()

        self.assertEqual( request3a.result(), 0 )
        self.assertEqual( request3b.result(), 0 )
        
        instance1.stop()
        instance2.stop()
        instance3.stop()

    def test21(self):
        """ test sum request, sockets """
        instance1 = ForTesting(self.exefile, channel_type="sockets")
        instance2 = ForTesting(self.exefile, channel_type="sockets")
        
        instance1.do_sleep(1, return_request=True)
        
        r1=instance1.echo_int(1, return_request=True)
        r2=instance2.echo_int(2, return_request=True)
        s=r1+r2
        r1=instance1.echo_int(2, return_request=True)
        r2=instance2.echo_int(3, return_request=True)
        m=r1*r2
        r1=instance1.echo_int(12, return_request=True)
        r2=instance2.echo_int(3, return_request=True)
        d=r1/r2
        
        self.assertEqual( s.result(), 3)
        self.assertEqual( m.result(), 6)
        self.assertEqual( d.result(), 4)
        
        instance1.stop()
        instance2.stop()

    def test21(self):
        """ some more tests of request expressions """
        instance1 = ForTesting(self.exefile)

        a=[10,30,15] | units.m
        b=[1,3,5] | units.kg
        
        instance1.do_sleep(1, return_request=True)
        request1=instance1.echo_2_int(a , b, return_request=True)
        
        request2=(3*request1[1]/(2.*request1[0])+(55. | units.kg/units.m))
        self.assertEqual( request2.result(), (3*b/(2.*a)+(55. | units.kg/units.m)) )
        
        instance1.stop()

    def test22(self):
        """ tests of unpack """
        instance1 = ForTesting(self.exefile)

        a=[10,30,15] | units.m
        b=[1,3,5] | units.kg
        
        #~ instance1.do_sleep(1, return_request=True)
        a_,b_=instance1.echo_2_int(a , b, return_request=True)
        
        self.assertEqual( (3*b_/(2.*a_)+(55. | units.kg/units.m)).result(), (3*b/(2.*a)+(55. | units.kg/units.m)) )
        
        instance1.stop()

    def test23(self):
        """ tests of unpack """
        instance1 = ForTestingInterface(self.exefile)

        a=[10,30,15]
        b=[1,3,5]
        
        #~ instance1.do_sleep(1, return_request=True)
        res=instance1.echo_2_int.asynchronous(a,b)
        #~ res=res['int_out1']
        a_,b_, err= res
      

        self.assertEqual( a,a_.result() )
        self.assertEqual( b,b_.result() )
        
        instance1.stop()

    def test24(self):
        """ more test of pool """
        instance1 = ForTesting(self.exefile)
        instance2 = ForTesting(self.exefile)
        
        r1=instance1.echo_int(1, return_request=True)
        r2=instance1.echo_int(2, return_request=True)
        r3=instance2.echo_int(3, return_request=True)
        r4=instance2.echo_int(4, return_request=True)
        
        p1=r1.join(r3)
        p2=r2.join(r4)
        
        p3=p1.join(p2)
        
        self.assertTrue(p3 is p1)
        
        p3.waitall()
        
        self.assertEqual(r1.result(), 1)
        self.assertEqual(r2.result(), 2)
        self.assertEqual(r3.result(), 3)
        self.assertEqual(r4.result(), 4)
        
        instance1.stop()
        instance2.stop()

    def test25(self):
        """ more test of pool: calls of same code """
        from amuse.rfi.async_request import AsyncRequestsPool
        instance1 = ForTesting(self.exefile)
        
        r1=instance1.do_sleep(1, return_request=True)
        r2=instance1.echo_int(2, return_request=True)
        
        p1=AsyncRequestsPool()
        r1.wait()
        r2.wait()
        p1.add_request(r1)
        p1.add_request(r2)
        
        #~ p1=r1.join(r2)
        
        p1.waitall()
        
        self.assertEqual(r2.result(), 2)
        
        instance1.stop()

    def test30(self):
        """ test a grid attribute request """
        instance1 = ForTesting(self.exefile, redirection="none")
        self.assertEqual(instance1.grid.x, numpy.arange(1,11) |units.m)
        instance1.do_sleep(1, return_request=True)
        t1=time.time()
        request=instance1.grid.request.x
        t2=time.time()
        self.assertLess(t2-t1, 0.5)
        self.assertEqual(request.result(), numpy.arange(1,11) | units.m)
        t2=time.time()
        self.assertGreater(t2-t1, 1.)
        
    def test31(self):
        """ test a grid attribute request, subgrids """
        instance1 = ForTesting(self.exefile, redirection="none")
        self.assertEqual(instance1.grid.x, numpy.arange(1,11) |units.m)
        instance1.do_sleep(1, return_request=True)
        t1=time.time()
        request=instance1.grid[:5].request.x
        request2=instance1.grid[5:].request.x
        t2=time.time()
        self.assertLess(t2-t1, 0.5)
        self.assertEqual(request.result(), numpy.arange(1,6) | units.m)
        self.assertEqual(request2.result(), numpy.arange(6,11) | units.m)
        t2=time.time()
        self.assertGreater(t2-t1, 1.)

    def test32(self):
        """ test a grid attribute request setter """
        instance1 = ForTesting(self.exefile, redirection="none")
        instance1.grid.x=(66.+numpy.arange(1,11)) |units.m
        self.assertEqual(instance1.grid.x, (66.+numpy.arange(1,11)) |units.m)

        t1=time.time()
        instance1.do_sleep(1, return_request=True)
        instance1.grid.request.x=(11.+numpy.arange(1,11)) |units.m
        t2=time.time()
        self.assertLess(t2-t1, 0.5)
        instance1.async_request.wait()
        t2=time.time()
        self.assertGreater(t2-t1, 1.)
        t1=time.time()
        self.assertEqual(instance1.grid.x, (11.+numpy.arange(1,11)) | units.m)
        t2=time.time()
        self.assertLess(t2-t1, 0.5)

    def test33(self):
        """ test a grid attribute request, subgrids """
        instance1 = ForTesting(self.exefile, redirection="none")
        self.assertEqual(instance1.grid.x, numpy.arange(1,11) |units.m)

        t1=time.time()
        instance1.do_sleep(1, return_request=True)
        instance1.grid[::2].request.x=(11.+numpy.arange(1,11,2)) |units.m
        t2=time.time()
        self.assertLess(t2-t1, 0.5)
        instance1.async_request.wait()
        t2=time.time()
        self.assertGreater(t2-t1, 1.)
        self.assertEqual(instance1.grid.x[::2], (11.+numpy.arange(1,11,2)) | units.m)
        self.assertEqual(instance1.grid.x[1::2], (numpy.arange(2,11,2)) | units.m)

    def test34(self):
        """ test a grid attribute request, subgrids """
        instance1 = ForTesting(self.exefile, redirection="none")
        grid=instance1.grid.copy()
        request=instance1.grid.request.x
        self.assertEqual(request.result(), numpy.arange(1,11) | units.m)

    def test35(self):
        """ test a grid attribute request setter with state"""
        instance1 = ForTestingWithState(self.exefile, redirection="none")
        t1=time.time()
        instance1.do_sleep(1, return_request=True)
        self.assertEqual(instance1.get_name_of_current_state(), '1')
        instance1.grid.request.x=(11.+numpy.arange(1,11)) |units.m
        self.assertEqual(instance1.get_name_of_current_state(), '2')
        t2=time.time()
        self.assertGreater(t2-t1, 1.)   # first time, state calls dummy (blocking) -> wait

        t1=time.time()
        instance1.do_sleep(1, return_request=True)
        instance1.grid.request.x=(12.+numpy.arange(1,11)) |units.m
        t2=time.time()
        self.assertLess(t2-t1, 0.5)  # second time should be less
                
        instance1.async_request.wait()
        t2=time.time()
        self.assertGreater(t2-t1, 1.)
        t1=time.time()
        self.assertEqual(instance1.grid.x, (12. +numpy.arange(1,11)) | units.m)
        t2=time.time()
        self.assertLess(t2-t1, 0.5)

    def test36(self):
        """ more state tests"""
        instance1 = ForTestingWithState(self.exefile, redirection="none")
        self.assertEqual(instance1.get_name_of_current_state(), '1')
        # this documents current behaviour:
        instance1.dummy(return_request=True)
        self.assertEqual(instance1.get_name_of_current_state(), '1')
        instance1.async_request.wait() 
        self.assertEqual(instance1.get_name_of_current_state(), '2')
        # ie state changes upon completion of call at wait. This is 
        # sort of ok, alternatively state could be changed immediately...


class TestASyncDistributed(TestASync):

    @classmethod
    def setup_class(cls):
        cls.check_not_in_mpiexec()
        super(TestASyncDistributed, cls).setup_class()
        cls.distinstance = cls.new_instance_of_an_optional_code(DistributedAmuse)#, redirection='none')
        cls.distinstance.parameters.debug = False

        #~ print "Resources:"
        #~ print cls.distinstance.resources

        pilot = Pilot()
        pilot.resource_name='local'
        pilot.node_count=1
        pilot.time= 2|units.hour
        pilot.slots_per_node=8
        pilot.label='local'
        cls.distinstance.pilots.add_pilot(pilot)
        #~ print "Pilots:"
        #~ print cls.distinstance.pilots

        #~ print "Waiting for pilots"
        cls.distinstance.wait_for_pilots()
        cls.distinstance.use_for_all_workers()

    @classmethod
    def tearDown(cls):
        #~ print "Stopping distributed code"
        cls.distinstance.stop()

    @classmethod
    def check_not_in_mpiexec(cls):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 
        
        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDI_CONTROL_FD' in os.environ:
            return
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            cls.skip('cannot run the socket tests under hydra process manager')


