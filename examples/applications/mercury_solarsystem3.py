import numpy
from amuse.community.mercury.interface import MercuryWayWard
from amuse.community.sse.interface import SSE
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
    timerange = units.day(numpy.arange(0.*365.25, 10000000 * 365.25, 92000))
    t_end = timerange[-1]
    gd = MercuryWayWard()
    gd.initialize_code()
    gd.central_particle.add_particles(sun)
    gd.orbiters.add_particles(planets)
    gd.commit_particles()

    se = SSE()
    #se.initialize_code()
    se.commit_parameters()
    se.particles.add_particles(sun)
    se.commit_particles()
    channelp = gd.orbiters.new_channel_to(planets)
    channels = se.particles.new_channel_to(sun)

    for time in timerange:
        if ((time.number % 920000) == 0):
            print time.value_in(units.yr)
            print sun.mass
        err = gd.evolve_model(time)
        channelp.copy()
        planets.savepoint(time)
        err = se.evolve_model(time + (12.32e9*365.25|units.day))
        channels.copy()
        gd.central_particle.mass = sun.mass


    gd.stop()
    se.stop()

    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        plot(x, y,'.')

    native_plot.show()

if __name__=='__main__':
    planetplot()
