import numpy
from amuse.community.mercury.interface import MercuryWayWard
from amuse.community.sse.interface import SSE
from amuse.ext.solarsystem import new_solar_system_for_mercury
from amuse.units import units
from amuse.units.quantities import VectorQuantity

from amuse.plot import *

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def planetplot():
    sun, planets = new_solar_system_for_mercury()

    initial = 12.2138 | units.Gyr
    final  =  12.3300 | units.Gyr
    step = 10000.0 | units.yr

    timerange = VectorQuantity.arange(initial, final, step)
    gd = MercuryWayWard()
    gd.initialize_code()
    #gd.stopping_conditions.timeout_detection.disable()
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
        err = gd.evolve_model(time-initial)
        channelp.copy()
        #planets.savepoint(time)
        err = se.evolve_model(time)
        channels.copy()
        gd.central_particle.mass = sun.mass
        print sun[0].mass.value_in(units.MSun), time.value_in(units.Myr), planets[4].x.value_in(units.AU),  planets[4].y.value_in(units.AU),planets[4].z.value_in(units.AU)

    gd.stop()
    se.stop()

    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        plot(x, y,'.')
        native_plot.gca().set_aspect('equal')

    native_plot.show()

if __name__ == "__main__":
    planetplot()
