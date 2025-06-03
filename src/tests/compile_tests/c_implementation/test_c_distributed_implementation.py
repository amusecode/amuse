import os
from amuse.rfi.core import *

from compile_tests.c_implementation.test_c_implementation import (
        TestCImplementationInterface)


class TestCDistributedImplementationInterface(TestCImplementationInterface):

    @classmethod
    def setup_class(cls):
        # Importing here so that the file can be imported by pytest even if
        # we don't have amuse.community.distributed available
        try:
            from amuse.community.distributed.interface import DistributedAmuse, Pilot
        except ImportError:
            cls.skip('Distributed AMUSE is not installed')

        cls.check_not_in_mpiexec()
        super(TestCDistributedImplementationInterface, cls).setup_class()
        # ~ print "Setting up distributed code"
        # instance = DistributedAmuse(redirection='none')
        cls.distinstance = cls.new_instance_of_an_optional_code(DistributedAmuse)  # , redirection='none')
        cls.distinstance.parameters.debug = False

        # ~ print "Resources:"
        # ~ print cls.distinstance.resources

        pilot = Pilot()
        pilot.resource_name = 'local'
        pilot.node_count = 1
        pilot.time = 2 | units.hour
        pilot.slots_per_node = 2
        pilot.label = 'local'
        cls.distinstance.pilots.add_pilot(pilot)
        # ~ print "Pilots:"
        # ~ print cls.distinstance.pilots

        # ~ print "Waiting for pilots"
        cls.distinstance.wait_for_pilots()
        cls.distinstance.use_for_all_workers()

    @classmethod
    def tearDown(cls):
        # ~ print "Stopping distributed code"
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
