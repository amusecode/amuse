"""
Example of usage of the distributed code based on the 'sunandearth' simple
example.
Evolves the dynamic evolution of the earth around the sun.
"""

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import quantities

from amuse.community.hermite.interface import Hermite

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

from matplotlib import pyplot

from amuse import datamodel

import webbrowser


def start_distributed_amuse():
    print "Creating distributed amuse"
    distributed_amuse = DistributedAmuse(redirection='none')
    distributed_amuse.parameters.debug = True
    distributed_amuse.parameters.webinterface_port = 4556

    distributed_amuse.use_for_all_workers()

    # open the address of the webinterface in a brower window
    webbrowser.open(distributed_amuse.get_webinterface_url())

    # Add some resources
    # resource = Resource()
    # resource.name='some.machine'
    # resource.location="user@fs0.das4.cs.vu.nl"
    # resource.scheduler_type="sge"
    # resource.amuse_dir="/home/user/amuse"
    # distributed_amuse.resources.add_resource(resource)

    print "Resources:"
    print distributed_amuse.resources

    # Claim nodes on the resources. In this example simply the "local" machine
    pilot = Pilot()
    pilot.resource_name = 'local'
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 22
    pilot.label = 'local'
    distributed_amuse.pilots.add_pilot(pilot)

    print "Pilots:"
    print distributed_amuse.pilots

    print "Waiting for pilots"
    distributed_amuse.wait_for_pilots()

    print "setting distributed as default channel"
    distributed_amuse.use_for_all_workers()

    return distributed_amuse


def new_system_of_sun_and_earth():
    stars = datamodel.Particles(2)

    sun = stars[0]
    sun.mass = 1.0 | units.MSun
    sun.position = (0.0, 0.0, 0.0) | units.m
    sun.velocity = (0.0, 0.0, 0.0) | (units.m / units.s)
    sun.radius = 1.0 | units.RSun

    earth = stars[1]
    earth.mass = 5.9736e24 | units.kg
    earth.radius = 6371.0 | units.km
    earth.position = (149.5e6, 0.0, 0.0) | units.km
    earth.velocity = (0.0, 29800, 0.0) | (units.m / units.s)

    return stars


def simulate_system_until(particles, end_time):
    convert_nbody = nbody_system.nbody_to_si(
        1.0 | units.MSun, 149.5e6 | units.km)

    instance = Hermite(convert_nbody)
    instance.parameters.epsilon_squared = 0.0 | units.AU**2
    instance.particles.add_particles(particles)

    t0 = 0 | units.day
    dt = 10 | units.day
    t = t0
    earth = instance.particles[1]

    x_values = quantities.AdaptingVectorQuantity()
    y_values = quantities.AdaptingVectorQuantity()

    while t < end_time:
        instance.evolve_model(t)

        x_values.append(earth.x)
        y_values.append(earth.y)

        t += dt

    instance.stop()

    return x_values, y_values


def plot_track(x, y):

    figure = pyplot.figure(figsize=(5, 5))
    plot = figure.add_subplot(1, 1, 1)

    x_values_in_AU = x.value_in(units.AU)
    y_values_in_AU = y.value_in(units.AU)

    plot.plot(x_values_in_AU, y_values_in_AU, color="b")

    plot.set_xlim(-1.5, 1.5)
    plot.set_ylim(-1.5, 1.5)

    plot.set_xlabel('x (AU)')
    plot.set_ylabel('y (AU)')

    pyplot.show()


if __name__ in ('__main__', '__plot__'):
    distributed_amuse = start_distributed_amuse()

    particles = new_system_of_sun_and_earth()
    x, y = simulate_system_until(particles, 20 | units.yr)
    plot_track(x, y)

    distributed_amuse.stop()
