from amuse.support.interface import InCodeComponentImplementation

from amusetest import TestWithMPI, TestCase

import numpy
import sys
import os
import time
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi import python_code
from amuse.rfi.core import *
from amuse.rfi.async_request import AsyncRequestsPool

from compile_tests import test_python_implementation
from compile_tests import test_python_implementation_mpi


class TestInterfaceSockets(test_python_implementation.TestInterface):
    def setUp(self):
        self.check_not_in_mpiexec()
        super(test_python_implementation.TestInterface, self).setUp()

    def check_not_in_mpiexec(self):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 

        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            self.skip('cannot run the socket tests under hydra process manager')

    def ForTesting(self, **options):
        options["worker_dir"] = self.get_path_to_results()
        options["channel_type"] = "sockets"
        return test_python_implementation.ForTesting(**options)

    def ForTestingInterface(self, **options):
        options["worker_dir"] = self.get_path_to_results()
        options["channel_type"] = "sockets"
        return test_python_implementation.ForTestingInterface(**options)

    def test2(self):
        implementation = test_python_implementation.ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, test_python_implementation.ForTestingInterface)

        input_message = python_code.SocketMessage(0, 10, 1)
        input_message.ints = [1]

        output_message = python_code.SocketMessage(0, 10, 1)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 1)
        self.assertEqual(len(output_message.doubles), 1)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(output_message.doubles[0], 0.0)

    def test3(self):
        implementation = test_python_implementation.ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, test_python_implementation.ForTestingInterface)

        input_message = python_code.SocketMessage(0, 11, 1)
        input_message.ints = [1]
        input_message.doubles = [12.0]

        output_message = python_code.SocketMessage(0, 10, 1)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 1)
        self.assertEqual(len(output_message.doubles), 0)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(implementation.masses[1], 12.0)

    def test4(self):
        implementation = test_python_implementation.ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, test_python_implementation.ForTestingInterface)

        input_message = python_code.SocketMessage(0, 11, 4)
        input_message.ints = [1, 2, 3, 4]
        input_message.doubles = [12.0, 13.0, 14.0, 15.0]

        output_message = python_code.SocketMessage(0, 10, 4)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 4)
        self.assertEqual(len(output_message.doubles), 0)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(output_message.ints[3], 0)
        self.assertEqual(implementation.masses[1], 12.0)
        self.assertEqual(implementation.masses[2], 13.0)
        self.assertEqual(implementation.masses[3], 14.0)
        self.assertEqual(implementation.masses[4], 15.0)

    def test8(self):
        implementation = test_python_implementation.ForTestingImplementation()
        x = python_code.PythonImplementation(implementation, test_python_implementation.ForTestingInterface)

        input_message = python_code.SocketMessage(0, 12, 1)
        input_message.ints = [20]

        output_message = python_code.SocketMessage(0, 10, 1)

        x.handle_message(input_message, output_message)

        self.assertEqual(len(output_message.ints), 2)
        self.assertEqual(output_message.ints[0], 0)
        self.assertEqual(output_message.ints[1], 20)

    def test26(self):
        # This fails with a SocketChannel because with a SocketChannel the workers
        # are started using mpiexec and connected to using TCP, and then this test
        # tries to open an MPI connection using the TCP connection as a side channel.
        # OpenMPI requires a separate daemon to be running for that, and we make no
        # effort to start it. Not sure why we're testing this in the first place, to
        # be investigated later - LV.
        self.skip("Skip because it fails with OpenMPI")

    def test27(self):
        self.skip("Only supported for mpi channel")


class TestInterfaceSocketsMPI(test_python_implementation_mpi.TestInterface):
    def setUp(self):
        self.check_not_in_mpiexec()
        super(test_python_implementation_mpi.TestInterface, self).setUp()

    def check_not_in_mpiexec(self):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 

        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            self.skip('cannot run the socket tests under hydra process manager')

    def ForTesting(self, **options):
        options["worker_dir"] = self.get_path_to_results()
        options["channel_type"] = "sockets"
        return test_python_implementation_mpi.ForTesting(**options)

    def ForTestingInterface(self, **options):
        options["worker_dir"] = self.get_path_to_results()
        options["channel_type"] = "sockets"
        return test_python_implementation_mpi.ForTestingInterface(**options)
