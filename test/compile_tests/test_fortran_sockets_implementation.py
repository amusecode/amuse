from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.test import compile_tools
import os
import time
import shlex

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi import channel
from amuse.rfi.core import *

from test_fortran_implementation import codestring, ForTestingInterface


class ForTesting(InCodeComponentImplementation):
    
    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)


class TestInterface(TestWithMPI):
        
    def setUp(self):
        super(TestInterface, self).setUp()
        print "building"
        self.check_can_compile_modules()
        self.check_fortran_version()
        self.check_not_in_mpiexec()
        self.exefile=compile_tools.build_fortran_worker(codestring, self.get_path_to_results(), ForTestingInterface)
        
    def check_fortran_version(self):
        pass
        
    def check_not_in_mpiexec(self):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 
        
        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDI_CONTROL_FD' in os.environ:
            return # can run in modern mpiexec.hydra
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            self.skip('cannot run the socket tests under hydra process manager')
    
        
    def test1(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        int_out, error = instance.echo_int(10)
        del instance
        self.assertEquals(int_out, 10)
        self.assertEquals(error, 0)
        
    def test2(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.echo_double(4.0)
        del instance
        self.assertEquals(out, 4.0)
        self.assertEquals(error, 0)
        
        
    def test3(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        input = [1,2,3,4]
        output, errors = instance.echo_int(input)
        del instance
        self.assertEquals(len(errors),4)
        for actual, expected in zip(output, input):
            self.assertEquals(actual, expected)
            
    def test4(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        input = [1.0,2.1,3.3,4.2]
        output, errors = instance.echo_double(input)
        del instance
        self.assertEquals(len(errors),4)
        for actual, expected in zip(output, input):
            self.assertEquals(actual, expected)
            
        
    def test5(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.echo_float(4.0)
        del instance
        self.assertEquals(out, 4.0)
        self.assertEquals(error, 0)
        
    def test6(self):
        
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.echo_string("abc")
        del instance
        self.assertEquals(error, 0)
        self.assertEquals(out, "abc")

    def test7(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.echo_string(["abc","def"])
        del instance
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(out[0], "abc")
        self.assertEquals(out[1], "def")

    def test8(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out1, out2, error = instance.echo_strings("abc","def")
        del instance
        
        self.assertEquals(error, 0)
        self.assertEquals(out1, "Abc")
        self.assertEquals(out2, "Bef")
      
    def test9(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        str1_out, str2_out, error = instance.echo_strings(["abc", "def"], ["ghi", "jkl"])
        del instance
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(str1_out[0], "Abc")
        self.assertEquals(str1_out[1], "Aef")
        self.assertEquals(str2_out[0], "Bhi")
        self.assertEquals(str2_out[1], "Bkl")
      
    def test10(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out = instance.return_string("qwerty")
        out = instance.return_string("abcdefghi")
        
        instance.stop()
        del instance
        
        self.assertEquals(out, "abcdefghi")
        
    def test11(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.hello_string()
        del instance
        
        self.assertEquals(out, "hello")
        
    def test12(self):
        
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.echo_string_fixed_len("abc")
        del instance
        self.assertEquals(error, 0)
        self.assertEquals(out, "abc")
        
        
    def test13(self):
        instance = ForTestingInterface(self.exefile, debugger="none", channel_type="sockets")
        (output_ints, error) = instance.echo_array_with_result([4,5,6])
        instance.stop()
        print output_ints, error
        self.assertEquals(output_ints[0], 4)
        self.assertEquals(output_ints[1], 5)
        self.assertEquals(output_ints[2], 6)
        
        self.assertEquals(error[0], -1)
        self.assertEquals(error[1], -1)
        self.assertEquals(error[2], -1)
        
    
    def test14(self):
        
        for x in range(4):
            instance = ForTestingInterface(self.exefile, channel_type="sockets")
            int_out, error = instance.echo_int(10)
            instance.stop()
            self.assertEquals(int_out, 10)
            self.assertEquals(error, 0)
            
    def test15(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        (output_ints, error) = instance.echo_inout_array_with_result([4,5,6])
        instance.stop()
        self.assertEquals(output_ints[0], 14)
        self.assertEquals(output_ints[1], 15)
        self.assertEquals(output_ints[2], 16)
        
        self.assertEquals(error[0], 11)
        self.assertEquals(error[1], 11)
        self.assertEquals(error[2], 11)
        
    def test16(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        (output1, error1) = instance.echo_logical(True)
        (output2, error2) = instance.echo_logical(False)
        instance.stop()
        self.assertEquals(error1, 0)
        self.assertEquals(error2, 0)
        self.assertTrue(output1)
        self.assertFalse(output2)

    def test16b(self):
        instance = ForTesting(self.exefile, channel_type="sockets")
        output = instance.echo_logical([True, True,False, True, False])
        self.assertEquals(output, [True, True, False, True, False])

    def test16c(self):
        instance = ForTesting(self.exefile, redirection="none")
        output = instance.echo_logical2([True, True,False, True, False]*1024)
        self.assertEquals(output, [True, True, False, True, False]*1024)
                
    def xtest20(self):
        #
        # TURNED OFF support for redirection,
        # by default output is redirected to /dev/null
        # if you need file, use the support from your mpi implementation
        #
        if os.path.exists("pout.000"):
            os.remove("pout.000")
        if os.path.exists("perr.000"):
            os.remove("perr.000")
        
        x = ForTestingInterface(self.exefile, redirect_stderr_file = 'perr', redirect_stdout_file = 'pout', redirection="file", channel_type="sockets")
        x.print_string("abc")
        x.print_error_string("exex")
        x.stop()
        
        time.sleep(0.2)
        
        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "abc")
        
        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "exex")
        
        x = ForTestingInterface(self.exefile, redirect_stderr_file = 'perr', redirect_stdout_file = 'pout', redirection="file", channel_type="sockets")
        x.print_string("def")
        x.print_error_string("exex")
        x.stop()
        
        time.sleep(0.2)
        
        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "abc\n def") 
        
        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000","r") as f:
            content = f.read()
        self.assertEquals(content.strip(), "exex\n exex")
        
    def test35(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        out, error = instance.echo_string(["abc","def"]*100000)
        del instance
        
        self.assertEquals(error[0], 0)
        self.assertEquals(error[1], 0)
        self.assertEquals(out[-2], "abc")
        self.assertEquals(out[-1], "def")
        
    def test36(self):
        instance = ForTestingInterface(self.exefile, channel_type="sockets")
        N=255
        out, error = instance.echo_string("a"*N)
        del instance
        
        self.assertEquals(error, 0)
        self.assertEquals(out, "a"*N)
