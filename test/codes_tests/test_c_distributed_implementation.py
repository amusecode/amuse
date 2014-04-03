from amuse.support.interface import InCodeComponentImplementation

from amuse.community.distributed.interface import DistributedAmuse, Pilot

from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions

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

class TestCDistributedImplementationInterface(test_c_implementation.TestCImplementationInterface):

    def setUp(self):
        super(TestCDistributedImplementationInterface, self).setUp()
        print "Setting up distributed code"
        #instance = DistributedAmuse(redirection='none')
        self.distinstance = self.new_instance_of_an_optional_code(DistributedAmuse)#, redirection='none')
        self.distinstance.parameters.debug = False

        print "Resources:"
        print self.distinstance.resources

        pilot = Pilot()
        pilot.resource_name='local'
        pilot.node_count=1
        pilot.time= 2|units.hour
        pilot.slots_per_node=2
        pilot.label='local'
        self.distinstance.pilots.add_pilot(pilot)
        print "Pilots:"
        print self.distinstance.pilots

        print "Waiting for pilots"
        self.distinstance.wait_for_pilots()
        self.distinstance.set_as_default()

    def tearDown(self):
        print "Stopping distributed code"
        self.distinstance.stop()

    def test22(self):
        self.skip("this test uses mpi internals, skip here")

