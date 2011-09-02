import numpy
from amuse.community.sse.interface import SSE
from amuse.support.data import core
from amuse.units import units
from amuse.ext import solarsystem
from amuse.plot import *

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def plottillagb():
    sse = SSE()
    sse.commit_parameters()
    sun, planets = solarsystem.new_solar_system_for_mercury()
    sse.particles.add_particles(sun)
    sse.commit_particles()
    channel = sse.particles.new_channel_to(sun)
    channel.copy()
    
    timerange = units.Myr(numpy.arange(11000, 15000,10))
    masses = []|units.MSun

    for time in timerange:
        sse.evolve_model(time)
        channel.copy()
        masses.append(sse.particles[0].mass)
        #print time.value_in(units.yr), sse.particles[0].mass.value_in(units.MSun)
        
    sse.stop()
    plot(timerange, masses,'.')
    native_plot.show()

if __name__ in ("__main__", "__plot__"):
    plottillagb()
