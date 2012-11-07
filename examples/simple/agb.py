"""
Evolves an the sun in it's AGB star phase using SSE
"""
    
import numpy
from matplotlib import pyplot

from amuse.community.sse.interface import SSE
from amuse.units import units
from amuse.ext import solarsystem
from amuse import plot

from amuse import datamodel

def plottillagb():
    sun = datamodel.Particle(
        mass = 1 | units.MSun,
        radius = 1 | units.RSun
    )
    
    sse = SSE()
    sse.particles.add_particle(sun)
    
    channel_from_se_to_memory = sse.particles.new_channel_to(sun.as_set())
    channel_from_se_to_memory.copy()
    
    masses = []|units.MSun

    timerange = numpy.arange(11500, 13500,10) | units.Myr
    for time in timerange:
        sse.evolve_model(time)
        channel_from_se_to_memory.copy()
        masses.append(sun.mass)
        print time.as_quantity_in(units.Myr), sun.mass.as_quantity_in(units.MSun)
        
    sse.stop()
    
    plot.plot(timerange, masses,'.')
    pyplot.show()

if __name__ in ("__main__", "__plot__"):
    plottillagb()
