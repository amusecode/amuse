#!/usr/bin/python
import sys
import webbrowser
from amuse.lab import *
from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Reservation, Reservations

#Example on how to use the distributed code. This example should run most (if not all) existing scripts.
#This example only uses local resources to run on. Add remote resources to get it to do something interesting.

if (len(sys.argv) < 2):
    print "usage: amuse.sh distributed_amuse.py existing_amuse_script.py"
    sys.exit(1)

print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
instance.initialize_code()

#open the address of the webinterface in a brower window
#webbrowser.open(instance.get_webinterface_url())

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
reservation.slots_per_node=22
reservation.node_label='local'
instance.reservations.add_reservation(reservation)

print "Reservations:"
print instance.reservations

print "Waiting for reservations"
instance.wait_for_reservations()

print "Running script"

script = sys.argv[1]

sys.argv = sys.argv[1:]

execfile(script)

print "script done, stopping distributed code"

instance.stop()


