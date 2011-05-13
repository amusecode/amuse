import numpy
from amuse.community.mercury.interface import MercuryInterface, MercuryWayWard
from amuse.ext.solarsystem import solarsystem
from amuse.support.units import units
from amuse.support.data import values
from amuse.plot import *

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class Phase(object):
    def __init__(self):
        self.x = values.AdaptingVectorQuantity()
        self.y = values.AdaptingVectorQuantity()
        self.z = values.AdaptingVectorQuantity()
        
planetsstore = []

def makestore():
    for i in range(9):
        instance = Phase()
        planetsstore.append(instance)

def store(particles):
    for i in range(9):
        planetsstore[i].x.append(particles[i].x)
        planetsstore[i].y.append(particles[i].y)    
        planetsstore[i].z.append(particles[i].z)

def planetplot():
    s = solarsystem()
    makestore()
    sun, planets = s.new_solarsystem()
    instance = MercuryWayWard()
    instance.initialize_code()
    instance.commit_parameters()
    instance.central_particle.add_particles(sun)
    instance.orbiters.add_particles(planets)
    channel=sun.new_channel_to(instance.central_particle)
    channel.copy()
    instance.commit_particles()

    t_end=365.25*20.0 | units.day
    time=0|units.day
    channels = instance.orbiters.new_channel_to(planets)
    while time < t_end:
        time = time + (8 | units.day)
        err=instance.evolve_model(time)
        channels.copy()
        store(planets)

    instance.stop()
    
    for i in range(9):
        plot(planetsstore[i].x, planetsstore[i].y)
    native_plot.show()

if __name__=='__main__':
    planetplot()
