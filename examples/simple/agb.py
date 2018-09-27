"""
Evolves an the sun in it's AGB star phase using SSE
"""
from __future__ import print_function

import numpy
from matplotlib import pyplot

from amuse.community.sse.interface import SSE
from amuse.units import units
# from amuse.ext import solarsystem

from amuse import datamodel


def plottillagb():
    sun = datamodel.Particle(
        mass=1 | units.MSun,
        radius=1 | units.RSun
    )

    sse = SSE()
    sse.particles.add_particle(sun)

    channel_from_se_to_memory = sse.particles.new_channel_to(sun.as_set())
    channel_from_se_to_memory.copy()

    masses = [] | units.MSun

    timerange = numpy.arange(11500, 13500, 10) | units.Myr
    for time in timerange:
        sse.evolve_model(time)
        channel_from_se_to_memory.copy()
        masses.append(sun.mass)
        print(time.as_quantity_in(units.Myr),
              sun.mass.as_quantity_in(units.MSun))

    sse.stop()

    figure = pyplot.figure(figsize=(6, 6))

    subplot = figure.add_subplot(1, 1, 1)
    subplot.plot(timerange.value_in(units.Gyr),
                 masses.value_in(units.MSun), '.')
    subplot.set_xlabel('t (Gyr)')
    subplot.set_ylabel('mass (MSun)')

    pyplot.show()


if __name__ == '__main__':
    plottillagb()
