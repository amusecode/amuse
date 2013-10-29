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
        super(TestCSocketsImplementationInterface, self).setUp()
	#set sockets channel as default channel
        options.GlobalOptions.instance().override_value_for_option("channel_type", "sockets")

    def tearDown(self):
	del options.GlobalOptions.instance().overriden_options["channel_type"]
	pass

    def test22(self):
	print 'this test uses mpi internals, skip here'

