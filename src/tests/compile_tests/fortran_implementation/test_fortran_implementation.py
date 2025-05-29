import unittest
from amuse.support.interface import InCodeComponentImplementation

from amuse.support.testing.amusetest import TestWithMPI
from compile_tests.fortran_implementation.interface import ForTestingInterface
import os
from pathlib import Path
import time

from amuse.units import units
from amuse.rfi.core import *


class ForTesting(InCodeComponentImplementation):

    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)

    """
    def define_methods(self, object):
        object.add_method(
            'echo_int',
            (units.m,)
            ,
            (
                units.m,
                object.ERROR_CODE,
            )
        )    
    """


class TestInterface(TestWithMPI):

    @classmethod
    def skip_if_fortran_does_not_support_mpi(cls):
        try:
            from amuse import config
            is_configured = hasattr(config, 'compilers') and hasattr(config.compilers, 'fc_iso_c_bindings')
        except ImportError:
            is_configured = False

        if is_configured and config.compilers.fc_iso_c_bindings:
            return
        else:
            cls.skip("cannot run test as fortran does not support iso c bindings")

    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'f_worker')

    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        del instance
        self.assertEqual(int_out, 10)
        self.assertEqual(error, 0)

    def test2(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_double(4.0)
        del instance
        self.assertEqual(out, 4.0)
        self.assertEqual(error, 0)

    def test3(self):
        instance = ForTestingInterface(self.exefile)
        input = [1, 2, 3, 4]
        output, errors = instance.echo_int(input)
        del instance
        self.assertEqual(len(errors), 4)
        for actual, expected in zip(output, input):
            self.assertEqual(actual, expected)

    def test4(self):
        instance = ForTestingInterface(self.exefile)
        input = [1.0, 2.1, 3.3, 4.2]
        output, errors = instance.echo_double(input)
        del instance
        self.assertEqual(len(errors), 4)
        for actual, expected in zip(output, input):
            self.assertEqual(actual, expected)

    def test5(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_float(4.0)
        del instance
        self.assertEqual(out, 4.0)
        self.assertEqual(error, 0)

    def test6(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string("abc")
        del instance
        self.assertEqual(error, 0)
        self.assertEqual(out, "abc")

    def test7(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string(["abc", "def"])
        del instance

        self.assertEqual(error[0], 0)
        self.assertEqual(error[1], 0)
        self.assertEqual(out[0], "abc")
        self.assertEqual(out[1], "def")

    def test8(self):
        instance = ForTestingInterface(self.exefile)
        out1, out2, error = instance.echo_strings("abc", "def")
        del instance

        self.assertEqual(error, 0)
        self.assertEqual(out1, "Abc")
        self.assertEqual(out2, "Bef")

    def test9(self):
        instance = ForTestingInterface(self.exefile)
        str1_out, str2_out, error = instance.echo_strings(["abc", "def"], ["ghi", "jkl"])
        del instance

        self.assertEqual(error[0], 0)
        self.assertEqual(error[1], 0)
        self.assertEqual(str1_out[0], "Abc")
        self.assertEqual(str1_out[1], "Aef")
        self.assertEqual(str2_out[0], "Bhi")
        self.assertEqual(str2_out[1], "Bkl")

    def test10(self):
        instance = ForTestingInterface(self.exefile)
        out = instance.return_string("abc")
        del instance

        self.assertEqual(out, "abc")

    def test11(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.hello_string()
        del instance

        self.assertEqual(out, "hello")

    def test12(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_fixed_len("abc")
        del instance
        self.assertEqual(error, 0)
        self.assertEqual(out, "abc")

    def test13(self):
        instance = ForTestingInterface(self.exefile, debugger="none")
        (output_ints, error) = instance.echo_array_with_result([4, 5, 6])
        instance.stop()
        print(output_ints, error)
        self.assertEqual(output_ints[0], 4)
        self.assertEqual(output_ints[1], 5)
        self.assertEqual(output_ints[2], 6)

        self.assertEqual(error[0], -1)
        self.assertEqual(error[1], -1)
        self.assertEqual(error[2], -1)

    def test14(self):

        for x in range(4):
            instance = ForTestingInterface(self.exefile)
            int_out, error = instance.echo_int(10)
            instance.stop()
            self.assertEqual(int_out, 10)
            self.assertEqual(error, 0)

    def test15(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints, error) = instance.echo_inout_array_with_result([4, 5, 6])
        instance.stop()
        self.assertEqual(output_ints[0], 14)
        self.assertEqual(output_ints[1], 15)
        self.assertEqual(output_ints[2], 16)

        self.assertEqual(error[0], 11)
        self.assertEqual(error[1], 11)
        self.assertEqual(error[2], 11)

    def test16(self):
        instance = ForTestingInterface(self.exefile)
        (output1, error1) = instance.echo_logical(True)
        (output2, error2) = instance.echo_logical(False)
        instance.stop()
        self.assertEqual(error1, 0)
        self.assertEqual(error2, 0)
        self.assertTrue(output1)
        self.assertFalse(output2)

    def test16b(self):
        instance = ForTesting(self.exefile)
        output = instance.echo_logical([True, True, False, True, False]*256)
        self.assertEqual(output, [True, True, False, True, False]*256)

    def test16c(self):
        instance = ForTesting(self.exefile, redirection="none")
        output = instance.echo_logical2([True, True, False, True, False]*1024)
        self.assertEqual(output, [True, True, False, True, False]*1024)

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

        x = ForTestingInterface(self.exefile, redirect_stderr_file='perr', redirect_stdout_file='pout', redirection="file")
        x.print_string("abc")
        x.print_error_string("exex")
        x.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "abc")

        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "exex")

        x = ForTestingInterface(self.exefile, redirect_stderr_file='perr', redirect_stdout_file='pout', redirection="file")
        x.print_string("def")
        x.print_error_string("exex")
        x.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "abc\n def")

        self.assertTrue(os.path.exists("perr.000"))
        with open("perr.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "exex\n exex")

    def test17(self):
        self.check_for_mpi()
        instance = ForTestingInterface(self.exefile)
        (output1, error1) = instance.internal__get_message_polling_interval()
        error2 = instance.internal__set_message_polling_interval(1234)
        (output2, error3) = instance.internal__get_message_polling_interval()
        instance.stop()
        self.assertEqual(error1, 0)
        self.assertEqual(output1, 0)
        self.assertEqual(error2, 0)
        self.assertEqual(error3, 0)
        self.assertEqual(output2, 1234)

    def test18(self):
        self.check_for_mpi()
        self.skip_if_fortran_does_not_support_mpi()

        instance = ForTestingInterface(self.exefile)
        t0 = time.time()
        (output1, error1) = instance.internal__get_message_polling_interval()
        t1 = time.time()
        error2 = instance.internal__set_message_polling_interval(500 * 1000)
        t2 = time.time()
        (output2, error3) = instance.internal__get_message_polling_interval()
        time.sleep(0.1)
        (output2, error3) = instance.internal__get_message_polling_interval()
        t3 = time.time()
        instance.stop()
        self.assertEqual(error1, 0)
        self.assertEqual(output1, 0)
        self.assertEqual(error2, 0)
        self.assertEqual(error3, 0)
        self.assertEqual(output2, 500 * 1000)
        # ~ print t1 - t0, t3 - t2
        # ~ self.assertTrue((t3 - t2) > 0.25)

    def test19(self):
        self.check_for_mpi()
        instance1 = ForTestingInterface(self.exefile)
        instance2 = ForTestingInterface(self.exefile)
        portname, error = instance1.internal__open_port()
        self.assertTrue(len(portname) > 0)
        request1 = instance1.internal__accept_on_port.asynchronous(portname)
        request2 = instance2.internal__connect_to_port.asynchronous(portname)
        request1.wait()
        request2.wait()
        port_id1, error1 = request1.result()
        port_id2, error2 = request2.result()
        self.assertTrue(port_id1 >= 0)
        self.assertTrue(port_id2 >= 0)
        self.assertEqual(error1, 0)
        self.assertEqual(error2, 0)

    @unittest.skip
    def test31(self):
        import time
        instance = ForTestingInterface(self.exefile)
        N = 5000
        t1 = time.time()
        for i in range(N):
            res, err = instance.echo_int([i])
        t2 = time.time()
        print("1 time:", t2-t1, (t2-t1)/N)
        instance.stop()

        instance = ForTesting(self.exefile)
        N = 5000
        t1 = time.time()
        for i in range(N):
            res = instance.echo_int([i] | units.m)
        t2 = time.time()
        print("2 time:", t2-t1, (t2-t1)/N)
        instance.stop()

    def test32(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.get_element_status(numpy.arange(10))
        del instance

        self.assertEqual(out, ["dry"]*10)

    def test33(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.get_element_status(numpy.arange(100))
        del instance

        self.assertEqual(out, ["dry"]*100)

    def test34(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string(["abc"]*14)
        del instance

        self.assertEqual(out, ["abc"]*14)
        self.assertEqual(error, [0]*14)

    def test35(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string(["abc", "def"]*100000)
        del instance

        self.assertEqual(error[0], 0)
        self.assertEqual(error[1], 0)
        self.assertEqual(out[-2], "abc")
        self.assertEqual(out[-1], "def")

    def test36(self):
        instance = ForTestingInterface(self.exefile)
        N = 255
        out, error = instance.echo_string("a"*N)
        del instance

        self.assertEqual(error, 0)
        self.assertEqual(out, "a"*N)
