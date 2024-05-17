from pathlib import Path
import unittest

from amusetest import TestWithMPI
from compile_tests.c_implementation.interface import ForTestingInterface
from compile_tests.c_implementation.test_c_implementation import ForTesting
from amuse.support import exceptions

import os
import time
from amuse.units import units
from amuse.rfi.core import *


class TestCImplementationInterface(TestWithMPI):

    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'c_worker')

    def test1(self):
        instance = ForTestingInterface(self.exefile)
        int_out, error = instance.echo_int(10)
        instance.stop()
        self.assertEqual(int_out, 10)
        self.assertEqual(error, 0)

    def test2(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_double(4.0)
        instance.stop()
        self.assertEqual(out, 4.0)
        self.assertEqual(error, 0)

    def test3(self):
        instance = ForTestingInterface(self.exefile)
        input = [1, 2, 3, 4]
        output, errors = instance.echo_int(input)
        instance.stop()
        self.assertEqual(len(errors), 4)
        for actual, expected in zip(output, input):
            self.assertEqual(actual, expected)

    def test4(self):
        instance = ForTestingInterface(self.exefile)
        input = [1.0, 2.1, 3.3, 4.2]
        output, errors = instance.echo_double(input)
        instance.stop()
        self.assertEqual(len(errors), 4)
        for actual, expected in zip(output, input):
            self.assertEqual(actual, expected)

    def test5(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_float(4.0)
        instance.stop()
        self.assertEqual(out, 4.0)
        self.assertEqual(error, 0)

    def test6(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string("abc")
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out, "abc")

    def test7(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string(["abc", "def"])
        instance.stop()

        self.assertEqual(error[0], 0)
        self.assertEqual(error[1], 0)
        self.assertEqual(out[0], "abc")
        self.assertEqual(out[1], "def")

    def test8(self):
        instance = ForTestingInterface(self.exefile)
        out1, out2, error = instance.echo_strings("abc", "def")
        instance.stop()

        self.assertEqual(error, 0)
        self.assertEqual(out1, "def")
        self.assertEqual(out2, "abc")

    def test9(self):
        instance = ForTestingInterface(self.exefile)
        str1_out, str2_out, error = instance.echo_strings(["abc", "def"], ["ghi", "jkl"])
        instance.stop()

        self.assertEqual(error[0], 0)
        self.assertEqual(error[1], 0)
        self.assertEqual(str1_out[0], "ghi")
        self.assertEqual(str1_out[1], "jkl")
        self.assertEqual(str2_out[0], "abc")
        self.assertEqual(str2_out[1], "def")

    def xtest10(self):
        """test for ticket #74, 'running out of os file descriptors'


        Note: this test takes a very long time, to enable it
        remove the 'X' infront of the test name, to disable it
        add an 'X'.
        Also note: to test this, you best set the ulimit
        to a low number (but not too low), for example
        ulimit -n 400
        """
        for x in range(400):
            instance = ForTestingInterface(self.exefile)
            out, error = instance.echo_float(4.0)
            if x % 100 == 0:
                print("x:", x)
            instance.stop()

    def test11(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints,) = instance.echo_array([4, 5, 6])
        instance.stop()
        print(output_ints)
        self.assertEqual(output_ints[0], 4)
        self.assertEqual(output_ints[1], 5)
        self.assertEqual(output_ints[2], 6)

    def test12(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints, error) = instance.echo_array_with_result([4, 5, 6])
        instance.stop()
        self.assertEqual(output_ints[0], 4)
        self.assertEqual(output_ints[1], 5)
        self.assertEqual(output_ints[2], 6)

        self.assertEqual(error[0], -1)
        self.assertEqual(error[1], -1)
        self.assertEqual(error[2], -1)

    def test13(self):
        instance = ForTesting(self.exefile)
        self.assertRaises(exceptions.AmuseException, instance.echo_int, [-1, -2] | units.m,
            expected_message=
                "Error when calling 'echo_int' of a 'ForTesting', errorcode is -1")
        instance.stop()

    def test14(self):
        instance = ForTesting(self.exefile)
        self.assertRaises(exceptions.CodeException, lambda: instance.echo_int())
        old_id = instance.legacy_interface.echo_int.specification.id
        instance.legacy_interface.echo_int.specification.id = -9
        self.assertRaises(exceptions.CodeException, lambda: instance.echo_int(1 | units.m))
        instance.legacy_interface.echo_int.specification.id = old_id
        instance.echo_int(1 | units.m)
        instance.stop()

    def test15(self):
        instance = ForTesting(self.exefile)
        output_ints1, output_ints2 = instance.echo_2_int([1, 2], [3, 4])
        output_ints3, output_ints4 = instance.echo_2_int([1, 2, 3])
        output_ints5, output_ints6 = instance.echo_2_int([5], [0])
        output_ints7, output_ints8 = instance.echo_2_int([5])
        instance.stop()
        self.assertEqual(output_ints1[1], 2)
        self.assertEqual(output_ints2[0], 3)
        self.assertEqual(output_ints2[1], 4)
        for i in range(3):
            self.assertEqual(output_ints3[i], i + 1)
            self.assertEqual(output_ints4[i], 1)
        self.assertEqual(output_ints5[0], 5)
        self.assertEqual(output_ints6[0], 0)
        self.assertEqual(output_ints7[0], 5)
        self.assertEqual(output_ints8[0], 1)

    def test16(self):
        instance = ForTesting(self.exefile)

        self.assertRaises(exceptions.AmuseException, lambda: instance.echo_int([]))
        instance.stop()

    def test17(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints, error) = instance.echo_inout_array_with_result([4, 5, 6])
        instance.stop()
        self.assertEqual(output_ints[0], 14)
        self.assertEqual(output_ints[1], 15)
        self.assertEqual(output_ints[2], 16)

        self.assertEqual(error[0], 11)
        self.assertEqual(error[1], 11)
        self.assertEqual(error[2], 11)

    def test18(self):
        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_logical([True, False, True])
        instance.stop()
        self.assertEqual(out, [True, False, True])
        self.assertEqual(error, 0)

    def test19(self):
        instance = ForTestingInterface(self.exefile)
        print(3935559000370003845)
        int_out, error = instance.echo_long_long_int(3935559000370003845)
        instance.stop()
        self.assertEqual(int_out, 3935559000370003845)
        self.assertEqual(error, 0)

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

        x = ForTesting(self.exefile, redirect_stderr_file='perr', redirect_stdout_file='pout', redirection="file")
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

        x = ForTesting(self.exefile, redirect_stderr_file='pout', redirect_stdout_file='pout', redirection="file")
        x.print_string("def")
        x.print_error_string("exex")
        x.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists("pout.000"))
        with open("pout.000", "r") as f:
            content = f.read()
        self.assertEqual(content.strip(), "abc\ndef\nexex")

    def test21(self):
        instance = ForTestingInterface(self.exefile)
        (output1, error1) = instance.internal__get_message_polling_interval()
        error2 = instance.internal__set_message_polling_interval(1234)
        (output2, error3) = instance.internal__get_message_polling_interval()
        instance.internal__set_message_polling_interval(0)
        instance.stop()
        self.assertEqual(error1, 0)
        self.assertEqual(output1, 0)
        self.assertEqual(error2, 0)
        self.assertEqual(error3, 0)
        self.assertEqual(output2, 1234)

    def test22(self):
        self.check_for_mpi()
        instance = ForTestingInterface(self.exefile)
        t0 = time.time()
        (output1, error1) = instance.internal__get_message_polling_interval()
        t1 = time.time()
        error2 = instance.internal__set_message_polling_interval(500 * 1000)
        t2 = time.time()
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

    def test23(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_int(1)
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out, "echo")

    def test24(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_int(1, "abc")
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out, "abc")

    def test25(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_int([1, 2])
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out[0], "echo")
        self.assertEqual(out[1], "echo")

    def test26(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_int([1, 2], ["abc", "def"])
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out[0], "abc")
        self.assertEqual(out[1], "def")

    def test27(self):

        instance = ForTestingInterface(self.exefile)
        out, error = instance.echo_string_int([1, 2], "abc")
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out[0], "abc")
        self.assertEqual(out[1], "abc")

    def test28(self):
        instance = ForTestingInterface(self.exefile)
        out1, out2, error = instance.echo_string_two(["one", "two"], "three")
        instance.stop()
        self.assertEqual(error, 0)
        self.assertEqual(out1[0], "one")
        self.assertEqual(out1[1], "two")
        self.assertEqual(out2[0], "three")
        self.assertEqual(out2[1], "three")

    def test29(self):
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

    def test30(self):
        from amuse.support.interface import ConvertArgumentsException
        instance = ForTesting(self.exefile)
        self.assertRaises(ConvertArgumentsException, instance.dummy_3_int, 2, 3, i=1, expected_message="got multiple values for argument 'i' of method dummy_3_int")
        instance.stop()

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
