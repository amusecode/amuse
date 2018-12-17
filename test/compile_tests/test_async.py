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
from amuse.rfi import channel
from amuse.rfi.core import *

import test_c_implementation

from amuse.test import compile_tools


codestring=test_c_implementation.codestring+"""
#include <unistd.h>

int do_sleep(int in) {
    sleep(in);
    return 0;
}

int return_error() {
    return -1;
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
        function.result_type = 'int32'
        return function 


class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)

class TestASync(TestWithMPI):

    def setUp(self):
        super(TestASync, self).setUp()
        print "building...",
        self.check_can_compile_modules()
        try:
            self.exefile = compile_tools.build_worker(codestring, self.get_path_to_results(), ForTestingInterface)
        except Exception as ex:
            print ex
            raise
        print "done"
        
    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)

    def test2(self):
        instance = ForTestingInterface(self.exefile)
        request = instance.echo_int.asynchronous(10)
        self.assertEqual(request, instance.async_request)
        request.wait()
        int_out,error=request.result()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        instance.stop()

    def test3(self):
        instance = ForTestingInterface(self.exefile)
        request1 = instance.do_sleep.asynchronous(1)
        request2 = instance.echo_int.asynchronous(10)
        self.assertEqual(request2, instance.async_request)
        request2.wait()
        int_out,error=request2.result()
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        instance.stop()

    def test4(self):
        instance = ForTesting(self.exefile)
        request1 = instance.do_sleep(1, return_request=True)
        request2 = instance.echo_int(10, return_request=True)
        self.assertEqual(request2, instance.async_request)
        instance.async_request.wait()
        int_out=request2.result()
        self.assertEquals(int_out, 10)
        instance.stop()

    def test5(self):
        instance = ForTesting(self.exefile)
        instance.do_sleep(1, return_request=True)
        requests=[]
        for x in range(10):
            requests.append(instance.echo_int(x, return_request=True))
        instance.async_request.wait()
        for i,x in enumerate(requests):
            self.assertEquals(x.result(), i)
        instance.stop()

    def test6(self):
        instance = ForTesting(self.exefile)
        requests=[]
        for x in range(10):
            requests.append(instance.echo_int(x, return_request=True))
        instance.async_request.wait()
        for i,x in enumerate(requests):
            self.assertEquals(x.result(), i)
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
            self.assertEquals(x[0].result(), x[1])
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
            self.assertEquals(x[0].result(), x[1])
        instance1.stop()
        instance2.stop()
        self.assertTrue(t2-t1 < 2.)

    def test9(self):
        instance = ForTesting(self.exefile)
        for x in range(10):
            instance.echo_int(x, return_request=True)
        results=instance.async_request.results
        self.assertEquals(results, range(10))
        instance.stop()

    def test10(self):
        instance = ForTesting(self.exefile)
        r1=instance.do_sleep(1, return_request=True)
        r2=instance.return_error( return_request=True)
        r3=instance.echo_int(10, return_request=True)
        try:
            results=instance.async_request.results
        except Exception as ex:
            print ex
        print r1.is_result_available()
        print r2.is_result_available()
        print r2.result()
        try:
            print r3.is_result_available()
        except Exception as ex:
            print ex

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
        request2=channel.DependentASyncRequest(request1, fac)

        request2.wait()
        
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
            print index
            results[index]=result
            return result
        
        request1.add_result_handler(safe_result,(1,))
        
        def fac():
          return instance2.echo_int(results[1], return_request=True)

        #~ request2=instance2.echo_int(??, async_factory=fac)
        request2=channel.DependentASyncRequest(request1, fac)

        request2.wait()
        
        self.assertEqual( request2.result(), 10)
        
        instance1.stop()
        instance2.stop()


