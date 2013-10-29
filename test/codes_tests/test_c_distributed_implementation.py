from amuse.support.interface import InCodeComponentImplementation

from amuse.community.distributed.interface import DistributedAmuse, Reservation

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
	self.distinstance = DistributedAmuse(redirection='none')
	self.distinstance.initialize_code()

	print "Resources:"
	print self.distinstance.resources

	reservation = Reservation()
	reservation.resource_name='local'
	reservation.node_count=1
	reservation.time= 2|units.hour
	reservation.slots_per_node=2
	reservation.node_label='local'
	self.distinstance.reservations.add_reservation(reservation)
	print "Reservations:"
	print self.distinstance.reservations

	print "Waiting for reservations"
	self.distinstance.wait_for_reservations()

    def tearDown(self):
	print "Stopping distributed code"
	self.distinstance.stop()

    def test22(self):
	print 'this test uses mpi internals, skip here'

