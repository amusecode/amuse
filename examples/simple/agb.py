from amuse.community.sse.interface import SSE

from amuse.support.data import core
from amuse.support.units import units
from amuse.ext.solarsystem import Solarsystem
import numpy

def plottillagb():
    sse = SSE()
    sse.commit_parameters()
    sun, planets = Solarsystem.new_solarsystem()
    sse.particles.add_particles(sun)
    sse.commit_particles()
    channel = sse.particles.new_channel_to(sun)
    channel.copy()
    
    timerange = units.Myr(numpy.arange(10, 15000,100))
    for time in timerange:
        print time
        sse.evolve_model(time)
        channel.copy()
        #print sun.mass
        print sse.particles[0].mass

    sse.stop()

if __name__ == "__main__":
    plottillagb()
