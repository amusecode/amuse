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

Planet_enum = range(9)

class Phase(object):
    def __init__(self):
        self.x = values.AdaptingVectorQuantity()
        self.y = values.AdaptingVectorQuantity()
        self.z = values.AdaptingVectorQuantity()
        
planetsstore = [Phase() for each in Planet_enum]

def store(particles):
    for i in range(9):
        planetsstore[i].x.append(particles[i].x)
        planetsstore[i].y.append(particles[i].y)    
        planetsstore[i].z.append(particles[i].z)

def planetplot():
    s = solarsystem()
    sun, planets = s.new_solarsystem()
    instance = MercuryWayWard()
    instance.initialize_code()
    instance.commit_parameters()
    instance.central_particle.add_particles(sun)
    channel=sun.new_channel_to(instance.central_particle)
    channel.copy()

    instance.orbiters.add_particles(planets)
    instance.commit_particles()
    print instance.central_particle.mass
    print instance.central_particle
    t_end=365.25*1.0 | units.day
    time=0|units.day
    #channels = instance.orbiters.new_channel_to(planets)
    while time < t_end:
        time = time + (8 | units.day)
        err=instance.evolve_model(time)
        #channels.copy()
        #planets.savepoint(time)
        print instance.orbiters[0].mass
        print instance.orbiters[2].position
        #print planets[0].x
        #store(planets)

    instance.stop()

    import pdb; pdb.set_trace()
    #for i in range(9):
    #    plot(planetsstore[i].x, planetsstore[i].y)

    x = (planets[0].get_timeline_of_attribute("x"))
    y = (planets[0].get_timeline_of_attribute("y"))
    print x
    print map(lambda (t,X) : units.AU(X.value_in(units.AU)), x)
    plot((map(lambda (t,X) : units.AU(X.value_in(units.AU)), x)), 
         (map(lambda (t,Y) : units.AU(Y.value_in(units.AU)), y)))
    native_plot.show()

if __name__=='__main__':
    planetplot()
