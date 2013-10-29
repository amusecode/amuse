#!/usr/bin/python
import nose
from amuse.lab import *

from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Reservation, Reservations

#Simple script to run nosetests using the distributed code. Should work for any existing amuse test.
#This example only runs the tests on the local machine.

print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
instance.initialize_code()

#Add some resources
#resource = Resource()
#resource.name='DAS4-VU'
#resource.location="user@fs0.das4.cs.vu.nl"
#resource.scheduler_type="sge"
#resource.amuse_dir="/home/user/amuse"
#instance.resources.add_resource(resource)
print "Resources:"
print instance.resources

#Claim nodes on the resources. In this example simply the "local" machine
reservation = Reservation()
reservation.resource_name='local'
reservation.node_count=1
reservation.time= 2|units.hour
reservation.slots_per_node=2
reservation.node_label='local'
instance.reservations.add_reservation(reservation)
print "Reservations:"
print instance.reservations

print "Waiting for reservations"
instance.wait_for_reservations()

print "Running tests"

nose.run()

print "all tests done, stopping distributed code"

instance.stop()
