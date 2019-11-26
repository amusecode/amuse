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

from . import test_c_implementation

class TestCDistributedImplementationInterface(test_c_implementation.TestCImplementationInterface):

    @classmethod
    def setup_class(cls):
        cls.check_not_in_mpiexec()
        super(TestCDistributedImplementationInterface, cls).setup_class()
        #~ print "Setting up distributed code"
        #instance = DistributedAmuse(redirection='none')
        cls.distinstance = cls.new_instance_of_an_optional_code(DistributedAmuse)#, redirection='none')
        cls.distinstance.parameters.debug = False

        #~ print "Resources:"
        #~ print cls.distinstance.resources

        pilot = Pilot()
        pilot.resource_name='local'
        pilot.node_count=1
        pilot.time= 2|units.hour
        pilot.slots_per_node=2
        pilot.label='local'
        cls.distinstance.pilots.add_pilot(pilot)
        #~ print "Pilots:"
        #~ print cls.distinstance.pilots

        #~ print "Waiting for pilots"
        cls.distinstance.wait_for_pilots()
        cls.distinstance.use_for_all_workers()

    @classmethod
    def tearDown(cls):
        #~ print "Stopping distributed code"
        cls.distinstance.stop()

    @classmethod
    def check_not_in_mpiexec(cls):
        """
        The tests will fork another process, if the test run
        is itself an mpi process, the tests may fail. 
        
        For the hydra process manager the tests will fail.
        So skip the tests if we detect hydra
        """
        if 'HYDI_CONTROL_FD' in os.environ:
            return
        if 'HYDRA_CONTROL_FD' in os.environ or 'PMI_FD' in os.environ:
            cls.skip('cannot run the socket tests under hydra process manager')

    def test22(self):
        self.skip("this test uses mpi internals, skip here")

