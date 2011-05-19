import numpy
from amuse.community.mercury.interface import MercuryWayWard
#from amuse.community.hermite0.interface import Hermite
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
    timerange = units.yr(numpy.arange(185000,210000,100))#numpy.arange(0.*365.25, 200 , 960))
    t_end = timerange[-1]
    gd = MercuryWayWard(debugger='xterm')
    #gd = Hermite()
    gd.initialize_code()
    gd.stopping_conditions.timeout_detection.disable()

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
        err = gd.evolve_model(time)
        print err, time, planets[4].x.value_in(units.AU),  planets[4].y.value_in(units.AU),planets[4].z.value_in(units.AU)
        channelp.copy()
        planets.savepoint(time)
        #print planets[4].x.value_in(units.AU),\
        #    planets[4].y.value_in(units.AU),\
        #    planets[4].z.value_in(units.AU)
        err = se.evolve_model(time + (12.32e9|units.yr))
        channels.copy()
        gd.central_particle.mass = sun.mass
        

    gd.stop()
    se.stop()

    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        plot(x, y,'.')
        native_plot.gca().set_aspect('equal')

    native_plot.show()

if __name__=='__main__':
    planetplot()
