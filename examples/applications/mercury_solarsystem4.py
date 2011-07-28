import numpy
from amuse.community.mercury.interface import MercuryWayWard
from amuse.community.sse.interface import SSE
from amuse.ext.solarsystem import new_solar_system_for_mercury
from amuse.support.units import units
from amuse.support.data.values import VectorQuantity

from amuse.plot import *

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class SSEWithMassEvolve(SSE):
    def __init__(self, **options):
        SSE.__init__(self, convert_nbody = None, **options)

    def evolve_mass(self, mass):
        current_mass = self.particles[0].mass
        current_time = self.particles[0].age

        rc = 60.0 | (units.Gyr / units.MSun)
        
        print "find"

        while current_mass > mass:
            timestep = rc * (current_mass - mass)
            current_time += timestep
            #print timestep, current_time

            self.evolve_model(current_time)
            current_mass = self.particles[0].mass
            
        print "done"
        return current_time, current_mass

def setup_codes(sun, planets):
    gd = MercuryWayWard()#debugger='xterm')
    gd.initialize_code()
    gd.central_particle.add_particles(sun)
    gd.orbiters.add_particles(planets)
    gd.commit_particles()
    
    se = SSEWithMassEvolve()
    se.commit_parameters()
    se.particles.add_particles(sun)
    se.commit_particles()
    return gd, se

def testsse():
    sse = SSEWithMassEvolve()
    sse.commit_parameters()
    sun, planets = new_solar_system_for_mercury()
    sse.particles.add_particles(sun)
    sse.commit_particles()
    channel = sse.particles.new_channel_to(sun)
    channel.copy()
    
    massrange = units.MSun(numpy.arange(1, 0.8 ,-0.001))
    masses = []|units.MSun
    timerange = [] | units.Myr

    for mass in massrange:
        sse.evolve_mass(mass)
        t,m = sse.evolve_mass(mass)
        timerange.append(t)
        channel.copy()
        masses.append(sse.particles[0].mass)
        
    sse.stop()
    plot(massrange, timerange,'.')
    native_plot.show()

def polyevolve():
    sun, planets = Solarsystem.new_solarsystem()

    gd, se = setup_codes(sun, planets)

    channelp = gd.orbiters.new_channel_to(planets)
    channels = se.particles.new_channel_to(sun)

    prev_mass = sun[0].mass

    massrange = units.MSun(numpy.arange(0.9, 0.89999 ,-1e-9))
    masses = []|units.MSun
    timerange = [] | units.Myr

    for i, mass in enumerate(massrange):
        time, dummymass = se.evolve_mass(mass)
        if i==0: initialtime = time
        print time,mass,"evolving gd"
        gdtime = time-initialtime
        print gdtime
        err = gd.evolve_model(gdtime)
        channelp.copy()
        planets.savepoint(time)
        channels.copy()
        gd.central_particle.mass = sun[0].mass
        print sun[0].mass
        print sun[0].mass.value_in(units.MSun), time.value_in(units.Myr), planets[4].x.value_in(units.AU),  planets[4].y.value_in(units.AU),planets[4].z.value_in(units.AU)

    gd.stop()
    se.stop()

    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        t, z = planet.get_timeline_of_attribute_as_vector("z")
        plot3(x, y, z,'.')
        native_plot.gca().set_aspect('equal')

    native_plot.show()

if __name__=='__main__':
    testsse()
    #polyevolve()
