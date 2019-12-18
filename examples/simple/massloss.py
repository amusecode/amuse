# -*- coding: ascii -*-
"""
Evolves the sun and earth where the sun will lose mass every 220th step.
"""
from __future__ import print_function
import numpy
from amuse.community.hermite.interface import Hermite
# from amuse.community.sse.interface import SSE
from amuse import datamodel
from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import VectorQuantity
from amuse.plot import (
    plot, native_plot)


def simulate_massloss(time):
    return units.MSun(
        0.5 * (1.0 + 1.0 / (1.0 + numpy.exp((time.value_in(time.unit) - 70.0) / 15.)))
    )


if __name__ == "__main__":

    convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)

    particles = datamodel.Particles(2)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.position = [0.0, 0.0, 0.0] | units.AU
    sun.velocity = [0.0, 0.0, 0.0] | units.AU / units.yr
    sun.radius = 1.0 | units.RSun

    earth = particles[1]
    earth.mass = 5.9736e24 | units.kg
    earth.radius = 6371.0 | units.km
    earth.position = [0.0, 1.0, 0.0] | units.AU
    earth.velocity = [2.0 * numpy.pi, -0.0001, 0.0] | units.AU / units.yr

    instance = Hermite(convert_nbody)
    instance.particles.add_particles(particles)

    channelp = instance.particles.new_channel_to(particles)

    start = 0 | units.yr
    end = 150 | units.yr
    step = 10 | units.day

    timerange = VectorQuantity.arange(start, end, step)

    masses = [] | units.MSun

    for i, time in enumerate(timerange):
        instance.evolve_model(time)
        channelp.copy()
        particles.savepoint(time)
        if (i % 220 == 0):
            instance.particles[0].mass = simulate_massloss(time)
        masses.append(instance.particles[0].mass)

    instance.stop()

    particle = particles[1]

    t, pos = particle.get_timeline_of_attribute_as_vector("position")
    distances = pos.lengths().as_quantity_in(units.AU)

    plot(timerange, distances, timerange, masses)

    native_plot.show()
