from amuse.legacy.smallN.muse_dynamics_mpi import *
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data.core import Particle

# MPI Debugging
from amuse.legacy.support import channel
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
#channel.MessageChannel.DEBUGGER = channel.MessageChannel.GDB

myunits = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)
MEarth = 5.9742e24 | units.kg

mult = SmallN(myunits)
particle = Particle(1)
particle.mass = 1.0 | nbody_system.mass 
particle.x = 0.0 | nbody_system.length
particle.y = 0.0 | nbody_system.length
particle.z = 0.0 | nbody_system.length
nbody_speed = nbody_system.length / nbody_system.time
particle.vx = 0.0 | nbody_speed 
particle.vy = 0.0 | nbody_speed 
particle.vz = 0.0 | nbody_speed 
mult.add_particle(particle)

particle = Particle(2)
particle.mass = 1.0 * MEarth 
particle.x = 1.0 | units.AU
particle.y = 0.0 | units.AU
particle.z = -10.0 | units.AU
particle.vx = 0.0 | units.AU / units.yr
particle.vy = 0.0 | units.AU / units.yr
particle.vz = 1.0 | units.AU / units.yr
mult.add_particle(particle)

for k in range(0,2):
    print mult.get_particle_by_index(k)
print mult.get_total_energy()
mult.report_multiples(level=1)

mult.evolve(verbose=True)

for k in range(0,2):
    print mult.get_particle_by_index(k)
print mult.get_total_energy()
mult.report_multiples(level=1)

