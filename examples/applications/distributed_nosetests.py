#!/usr/bin/python
import nose


from amuse.units import units

from amuse.community.distributed.interface import (
        # DistributedAmuseInterface,
        DistributedAmuse,
        )
from amuse.community.distributed.interface import (
        # Resource,
        # Resources,
        Pilot,
        # Pilots,
        )

# Simple script to run nosetests using the distributed code. Should work for
# any existing amuse test.
# This example only runs the tests on the local machine.

print("Setting up distributed code")
instance = DistributedAmuse(redirection='none')
# instance.parameters.debug = True
# instance.parameters.webinterface_port = 4556
instance.commit_parameters()

# Add some resources
# resource = Resource()
# resource.name='DAS4-VU'
# resource.location="user@fs0.das4.cs.vu.nl"
# resource.scheduler_type="sge"
# resource.amuse_dir="/home/user/amuse"
# instance.resources.add_resource(resource)
print("Resources:")
print(instance.resources)

# Claim nodes on the resources. In this example simply the "local" machine
pilot = Pilot()
pilot.resource_name = 'local'
pilot.node_count = 1
pilot.time = 2 | units.hour
pilot.slots_per_node = 32
pilot.label = 'local'
instance.pilots.add_pilot(pilot)
print("Pilots:")
print(instance.pilots)

print("Waiting for pilots")
instance.wait_for_pilots()

print("setting distributed as default channel")
instance.use_for_all_workers()

print("Running tests")

nose.run()

print("all tests done, stopping distributed code")

instance.stop()
