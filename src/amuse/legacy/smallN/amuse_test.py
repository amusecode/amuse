from amuse.legacy.smallN.muse_dynamics_mpi import *
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data.core import Particle
from time import time
import math, sys

# MPI Debugging
from amuse.support.legacy import channel
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.GDB
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD

if __name__ == "__main__":
    test_list = [1, 2]

    myunits = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)
    MEarth = 5.9742e24 | units.kg

    mult = SmallN(myunits)
    
    if 1 in test_list:
        print "=== Test 1: An Earth-mass body moving past the Sun at low speed ==="
        particle = Particle(key=1)
        particle.mass = 1.0 | nbody_system.mass 
        particle.radius = 0.0 | nbody_system.length
        particle.x = 0.0 | nbody_system.length
        particle.y = 0.0 | nbody_system.length
        particle.z = 0.0 | nbody_system.length
        particle.vx = 0.0 | nbody_system.speed 
        particle.vy = 0.0 | nbody_system.speed 
        particle.vz = 0.0 | nbody_system.speed 
        mult.particles.add_particle(particle)

        particle = Particle(key=2)
        particle.mass = 1.0 * MEarth 
        particle.radius = 0.0 | nbody_system.length
        particle.x = 1.0 | units.AU
        particle.y = 0.0 | units.AU
        particle.z = -10.0 | units.AU
        particle.vx = 0.0 | units.AU / units.yr
        particle.vy = 0.0 | units.AU / units.yr
        particle.vz = 1.0 | units.AU / units.yr
        mult.particles.add_particle(particle)

        for k in range(0,2):
            print "==> Dump of Particle %d" % (k+1)
            print mult.particles[k]
        print "==> Total Energy is " , mult.total_energy.as_quantity_in(units.J)
        mult.report_multiples(level=1)

        print "==> Particles set up.  Evolving until a stable regime is reached."
        sys.stdout.flush()
        start_time = time()
        mult.evolve(verbose=True)
        end_time = time()
        sys.stdout.flush()
        print "==> Evolution complete.  Real time used: %.3f milliseconds." % ((end_time-start_time)*1000)
        print "==> NB: Evolve method is not necessarily synchronised to your output stream."
        sys.stdout.flush()

        for k in range(0,2):
            print "==> Dump of Particle %d" % (k+1)
            print mult.particles[k]
        print "==> Total Energy is " , mult.total_energy.as_quantity_in(units.J)
        mult.report_multiples(level=1)
        print "==> Interaction took %.2f years" % mult.get_time().value_in(units.yr)
        time_in_nbody = myunits.to_nbody(mult.get_time())
        print "==>   = %.2e nbody times" % time_in_nbody.number

        mult.reset_close_encounter()
        print "===================END TEST========================"
        sys.stdout.flush()


    if 2 in test_list:
        print "=== Test 2: An Earth-mass body orbiting the Sun ==="
        particle = Particle(key=1)
        particle.mass = 1.0 | nbody_system.mass 
        particle.radius = 0.0 | nbody_system.length
        particle.x = 0.0 | nbody_system.length
        particle.y = 0.0 | nbody_system.length
        particle.z = 0.0 | nbody_system.length
        particle.vx = 0.0 | nbody_system.speed
        particle.vy = 0.0 | nbody_system.speed 
        particle.vz = 0.0 | nbody_system.speed 
        mult.particles.add_particle(particle)

        particle = Particle(key=2)
        particle.mass = 1.0 * MEarth 
        particle.radius = 0.0 | nbody_system.length
        particle.x = 1.0 | units.AU
        particle.y = 0.0 | units.AU
        particle.z = 0.0 | units.AU
        particle.vx = 0.0 | units.AU / units.yr
        particle.vy = 2.0*math.pi | units.AU / units.yr
        particle.vz = 0.0 | units.AU / units.yr
        mult.particles.add_particle(particle)
    

        for k in range(0,2):
            print "==> Dump of Particle %d" % (k+1)
            print mult.particles[k]
        print "==> Total Energy is " , mult.total_energy.as_quantity_in(units.J)
        print "==> Total Energy is " , mult.get_total_energy()
        mult.report_multiples(level=1)

        print "==> Particles set up.  Evolving until a stable regime is reached."
        sys.stdout.flush()
        start_time = time()
        mult.evolve(verbose=True)
        end_time = time()
        sys.stdout.flush()
        print "==> Evolution complete.  Real time used: %.3f milliseconds." % ((end_time-start_time)*1000)
        print "==> NB: Evolve method is not necessarily synchronised to your output stream."
        sys.stdout.flush()

        for k in range(0,2):
            print "==> Dump of Particle %d" % (k+1)
            print mult.particles[k]
        print "==> Total Energy is " , mult.total_energy.as_quantity_in(units.J)
        mult.report_multiples(level=1)
        print "==> Interaction took %.2f years" % mult.get_time().value_in(units.yr)
        time_in_nbody = myunits.to_nbody(mult.get_time())
        print "==>   = %.2e nbody times" % time_in_nbody.number

        mult.reset_close_encounter()
        print "===================END TEST========================"
        sys.stdout.flush()


