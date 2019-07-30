from amuse.support.interface import InCodeComponentImplementation

from amuse.community.distributed.interface import DistributedAmuse, Pilot

from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions

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
        self.check_not_in_mpiexec()
        super(TestCDistributedImplementationInterface, self).setUp()
        #~ print "Setting up distributed code"
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
    def test22(self):
        self.skip("this test uses mpi internals, skip here")

