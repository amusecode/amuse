from amuse.legacy.smallN.muse_dynamics_mpi import *
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data.core import Particle
import math

# MPI Debugging
from amuse.legacy.support import channel
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.GDB
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD

if __name__ == "__main__":
    myunits = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)
    myunits.set_as_default()
    MEarth = 5.9742e24 | units.kg

    mult = SmallN(myunits)
    
    if True:
        particle = Particle(key=1)
        particle.mass = 1.0 | nbody_system.mass 
        particle.x = 0.0 | nbody_system.length
        particle.y = 0.0 | nbody_system.length
        particle.z = 0.0 | nbody_system.length
        nbody_speed = nbody_system.length / nbody_system.time
        particle.vx = 0.0 | nbody_speed 
        particle.vy = 0.0 | nbody_speed 
        particle.vz = 0.0 | nbody_speed 
        mult.add_particle(particle)

        particle = Particle(key=2)
        particle.mass = 1.0 * MEarth 
        particle.x = 1.0 | units.AU
        particle.y = 0.0 | units.AU
        particle.z = -10.0 | units.AU
        particle.vx = 0.0 | units.AU / units.yr
        particle.vy = 0.0 | units.AU / units.yr
        particle.vz = 1.0 | units.AU / units.yr
        mult.add_particle(particle)


#    p3 = Particle(key=3)
#    p3.mass = 1.0 | units.MSun
#    p4 = Particle(key=4)
#    p4.mass = 1.0 | units.MSun
#    binary = BinaryStar(key=11, component_particle1=p3, component_particle2=p4, \
#                    period = 1.0 | units.yr, eccentricity = 0.1)
#    mult.report_multiples(level=1)

        for k in range(0,2):
            print mult.get_particle_by_index(k)
        print mult.get_total_energy()
        mult.report_multiples(level=1)

        mult.evolve(verbose=True)

        for k in range(0,2):
            print mult.get_particle_by_index(k)
        print mult.get_total_energy()
        mult.report_multiples(level=1)
        print "Interaction took %.2f years" % mult.get_time().value_in(units.yr)
        time_in_nbody = myunits.to_nbody(mult.get_time())
        print "  = %.2e nbody times" % time_in_nbody.number

        mult.reset_close_encounter()

    particle = Particle(key=1)
    particle.mass = 1.0 | nbody_system.mass 
    particle.x = 0.0 | nbody_system.length
    particle.y = 0.0 | nbody_system.length
    particle.z = 0.0 | nbody_system.length
    nbody_speed = nbody_system.length / nbody_system.time
    particle.vx = 0.0 | nbody_speed 
    particle.vy = 0.0 | nbody_speed 
    particle.vz = 0.0 | nbody_speed 
    mult.add_particle(particle)

    particle = Particle(key=2)
    particle.mass = 1.0 * MEarth 
    particle.x = 1.0 | units.AU
    particle.y = 0.0 | units.AU
    particle.z = 0.0 | units.AU
    particle.vx = 0.0 | units.AU / units.yr
    particle.vy = 2.0*math.pi | units.AU / units.yr
    particle.vz = 0.0 | units.AU / units.yr
    mult.add_particle(particle)

    start_time = 100 | units.yr
    mult.set_time(start_time)

    for k in range(0,2):
        print mult.get_particle_by_index(k)
    print mult.get_total_energy()
    mult.report_multiples(level=1)

    #mult.evolve(super_verbose=True)
    mult.evolve()

    for k in range(0,2):
        print mult.get_particle_by_index(k)
    print mult.get_total_energy()
    mult.report_multiples(level=1)

    end_time = mult.get_time()
    duration = end_time - start_time
    print "Interaction took %.2f years" % duration.value_in(units.yr)
    time_in_nbody = myunits.to_nbody(duration)
    print "  = %.2e nbody times" % time_in_nbody.number
    mult.reset_close_encounter()


