from amuse.support.interface import InCodeComponentImplementation
from amusetest import TestWithMPI
from amuse.support import exceptions

from compile_tests.java_implementation.interface import ForTestingInterface

import os
import time
from amuse.units import units
from amuse.rfi.core import *


class ForTesting(InCodeComponentImplementation):

    def __init__(self, exefile, **options):
        InCodeComponentImplementation.__init__(self, ForTestingInterface(exefile, **options), **options)

    def define_methods(self, object):
        object.add_method(
            'echo_int',
            (units.m,),
            (
            units.m,
                object.ERROR_CODE,
            )
        )


class TestInterface(TestWithMPI):

    @classmethod
    def check_not_in_mpiexec(cls):
        """
        The tests will fork another process, if the test run 
        is itself an mpi process, the tests may fail.

        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDI_CONTROL_FD' in os.environ:
            return  # can run in modern mpiexec.hydra
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            cls.skip('cannot run the socket tests under mpi process manager')

    @classmethod
    def get_worker(cls):
        """
        Check if the worker is available. If it's not there,
        then we don't have Java available and skip the test.
        If it is there, return its location.
        """
        code_dir = os.path.dirname(__file__)
        worker_path = os.path.join(code_dir, "java_worker")
        if not os.path.exists(worker_path):
            cls.skip("java not enabled")
        return worker_path

    @classmethod
    def setup_class(cls):
        cls.check_not_in_mpiexec()
        cls.get_worker()
        cls.exefile = cls.get_worker()

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
        output, errors = instance.echo_int_array(input)
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
            instance.stop()

    def test11(self):
        instance = ForTestingInterface(self.exefile)
        (output_ints,) = instance.echo_array([4, 5, 6])
        instance.stop()
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
            expected_message="Error when calling 'echo_int' of a 'ForTesting', errorcode is -1")
        instance.stop()

    def test14(self):
        instance = ForTesting(self.exefile)
        self.assertRaises(exceptions.CodeException, lambda: instance.echo_int())
        instance.legacy_interface.echo_int.specification.id = -9
        self.assertRaises(exceptions.CodeException, lambda: instance.echo_int(1 | units.m))
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

        # self.assertRaises(exceptions.AmuseException, lambda : instance.echo_int([]))
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
        int_out, error = instance.echo_long_long_int(3935559000370003845)
        instance.stop()
        self.assertEqual(int_out, 3935559000370003845)
        self.assertEqual(error, 0)

    def test20(self):
        path = os.path.abspath(self.get_path_to_results())
        output = os.path.join(path, "output.txt")
        error = os.path.join(path, "error.txt")

        if os.path.exists(output):
            os.remove(output)
        if os.path.exists(error):
            os.remove(error)

        instance = ForTesting(self.exefile, redirect_stderr_file=error, redirect_stdout_file=output, redirection="file")
        instance.print_string("test_string_123")
        instance.print_error_string("error_string_123")
        instance.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists(output))
        print(f'Output at {path}')
        with open(output, "r") as f:
            content = f.read()
        self.assertTrue("test_string_123" in content.strip())

        self.assertTrue(os.path.exists(error))
        with open(error, "r") as f:
            content = f.read()
        # some times java generates "Picked up _JAVA_OPTIONS" message, so only test:
        self.assertTrue(content.strip().endswith("error_string_123"))

    def test21(self):
        path = os.path.abspath(self.get_path_to_results())
        output = os.path.join(path, "output.txt")
        error = os.path.join(path, "error.txt")

        if os.path.exists(output):
            os.remove(output)
        if os.path.exists(error):
            os.remove(error)

        instance = ForTesting(self.exefile, redirect_stderr_file=output, redirect_stdout_file=output, redirection="file")

        instance.print_string("abcdef")
        instance.print_error_string("&Hfecd")
        instance.stop()

        time.sleep(0.2)

        self.assertTrue(os.path.exists(output))
        with open(output, "r") as f:
            content = f.read()
        # some times java generates "Picked up _JAVA_OPTIONS" message, so only test:
        self.assertTrue("abcdef" in content)
        self.assertTrue("&Hfecd" in content)
