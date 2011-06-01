import numpy
from amuse.community.hermite0.interface import Hermite
from amuse.community.sse.interface import SSE
from amuse.ext.solarsystem import Solarsystem
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data.values import VectorQuantity
from amuse.support.data import core

from amuse.plot import *

def mass(time):
    return units.MSun(0.5*(1.0+1.0/(1.0+numpy.exp((time.value_in(time.unit)-70.0)/15.))))

def distance(x,y,z):
    return (x**2+y**2+z**2)**0.5

if __name__ == '__main__':

    convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
    
    particles = core.Particles(2)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.position = [0.0, 0.0, 0.0] | units.AU
    sun.velocity = [0.0, 0.0, 0.0] | units.AU / units.yr
    sun.radius = 1.0 | units.RSun
    
    earth = particles[1]
    earth.mass = 5.9736e24 | units.kg
    earth.radius = 6371.0 | units.km
    earth.position = [0.0, 1.0, 0.0] | units.AU
    earth.velocity = [2.0*numpy.pi, -0.0001, 0.0] | units.AU / units.yr
    
    instance = Hermite(convert_nbody)
    instance.initialize_code()
    instance.particles.add_particles(particles)
    instance.commit_particles()

    channelp = instance.particles.new_channel_to(particles)
    
    start = 0 |units.yr
    end = 150 | units.yr
    step = 10|units.day

    timerange = VectorQuantity.arange(start, end, step)

    masses = []|units.MSun

    for i, time in enumerate(timerange):
        instance.evolve_model(time)
        channelp.copy()
        particles.savepoint(time)
        if (i % 220 == 0):
            instance.particles[0].mass = mass(time)
        masses.append(instance.particles[0].mass)
 
    instance.stop()

    particle = particles[1]
    t, x = particle.get_timeline_of_attribute_as_vector("x")
    t, y = particle.get_timeline_of_attribute_as_vector("y")
    t, z = particle.get_timeline_of_attribute_as_vector("z")
    
    #plot3(x,y,z)
    plot(timerange,  distance(x,y,z), timerange, masses)

    #native_plot.gca().set_aspect('equal')
    
    native_plot.show()
