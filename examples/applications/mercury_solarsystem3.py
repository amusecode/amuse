import numpy
from amuse.community.mercury.interface import MercuryWayWard
from amuse.community.mesa.interface import MESA
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
    timerange = units.day(numpy.arange(0, 50000 * 365.25, 50))
    t_end = timerange[-1]
    gd = MercuryWayWard()
    gd.initialize_code()
    gd.central_particle.add_particles(sun)
    gd.orbiters.add_particles(planets)
    gd.commit_particles()

    se = MESA()
    se.initialize_code()
    se.parameters.RGB_wind_scheme = 1
    se.parameters.reimers_wind_efficiency = 1.0e6 # ridiculous, but instructive
    se.particles.add_particles(sun)

    channelp = gd.orbiters.new_channel_to(planets)
    channels = se.particles.new_channel_to(sun)

    for time in timerange:
        err = gd.evolve_model(time)
        channelp.copy()
        planets.savepoint(time)
        
        if (time.number % 100000==0):
            err = se.evolve_model(time+(5e7|units.day))
            channels.copy()
            gd.central_particle.mass = sun.mass
            sun.savepoint(time)
            print("\r {0} %%done".format(time/t_end * 100.0))
    gd.stop()

    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        plot(x, y,'.')

    native_plot.show()

if __name__=='__main__':
    planetplot()
