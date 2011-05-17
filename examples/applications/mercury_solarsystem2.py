import numpy
from amuse.community.mercury.interface import MercuryWayWard
from amuse.ext.solarsystem import Solarsystem
from amuse.support.units import units
from amuse.plot import *

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def planetplot():
    sun, planets = Solarsystem.new_solarsystem()
    timerange = units.day(numpy.arange(0, 120 * 365.25, 12))

    instance = MercuryWayWard()
    instance.initialize_code()
    instance.central_particle.add_particles(sun)
    instance.orbiters.add_particles(planets)
    instance.commit_particles()

    channels = instance.orbiters.new_channel_to(planets)

    for time in timerange:
        err = instance.evolve_model(time)
        channels.copy()
        planets.savepoint(time)

    instance.stop()

    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        plot(x, y,'.')

    native_plot.show()

if __name__=='__main__':
    planetplot()
