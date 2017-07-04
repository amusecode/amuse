from amuse.support.interface import InCodeComponentImplementation

from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions
from amuse.support import options

import subprocess
import os
import time
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi.tools import create_c
from amuse.rfi import channel
from amuse.rfi.core import *

import test_c_implementation

class TestCSocketsImplementationInterface(test_c_implementation.TestCImplementationInterface):

    def setUp(self):
        self.check_not_in_mpiexec()
        super(TestCSocketsImplementationInterface, self).setUp()
        #set sockets channel as default channel
        options.GlobalOptions.instance().override_value_for_option("channel_type", "sockets")

    def tearDown(self):
        del options.GlobalOptions.instance().overriden_options["channel_type"]

    def test22(self):
        self.skip("this test uses mpi internals, skip here")
        
    def test29(self):
        self.skip("this test uses mpi internals, skip here")
                     
    def check_not_in_mpiexec(self):
        """
        The tests will fork another process, if the test run 
        is itself an mpi process, the tests may fail.
                 
        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDI_CONTROL_FD' in os.environ:
            return # for now assume HYDI_CONTROL_FD is newer, and sockets will work!
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            self.skip('cannot run the socket tests under mpi process manager')
         
