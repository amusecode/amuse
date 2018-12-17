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

from amuse.community.distributed.interface import DistributedAmuse, Pilot

codestring=test_c_implementation.codestring+"""
#include <unistd.h>

int do_sleep(int in) {
    sleep(in);
    return 0;
}

int return_error(int * out) {
    *out=123;
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
        function.addParameter('out', dtype='int32', direction=function.OUT)
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
        r3=instance.echo_int(1, return_request=True)
        r4=instance.echo_int(2, return_request=True)
        self.assertRaises(Exception, instance.async_request.waitall,
            expected_message="Error when calling 'return_error' of a 'ForTesting', errorcode is -1")
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

class TestASyncDistributed(TestASync):

    def setUp(self):
        self.check_not_in_mpiexec()
        super(TestASyncDistributed, self).setUp()
        #instance = DistributedAmuse(redirection='none')
        self.distinstance = self.new_instance_of_an_optional_code(DistributedAmuse)#, redirection='none')
        self.distinstance.parameters.debug = False

        #~ print "Resources:"
        #~ print self.distinstance.resources

        pilot = Pilot()
        pilot.resource_name='local'
        pilot.node_count=1
        pilot.time= 2|units.hour
        pilot.slots_per_node=2
        pilot.label='local'
        self.distinstance.pilots.add_pilot(pilot)
        #~ print "Pilots:"
        #~ print self.distinstance.pilots

        #~ print "Waiting for pilots"
        self.distinstance.wait_for_pilots()
        self.distinstance.use_for_all_workers()

    def tearDown(self):
        #~ print "Stopping distributed code"
        self.distinstance.stop()

    def check_not_in_mpiexec(self):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 
        
        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDI_CONTROL_FD' in os.environ:
            return
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            self.skip('cannot run the socket tests under hydra process manager')


