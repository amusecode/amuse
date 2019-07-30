from amuse.lab import *
from interface import symple

n = 10

total_mass = n|units.MSun
length = 1|units.parsec
converter = nbody_system.nbody_to_si(total_mass, length)

# Select a gravity code.

g = symple(convert_nbody=converter, redirection='none')
g.initialize_code()
g.parameters.set_defaults()

# Make an N-body system.

p = new_plummer_model(n, convert_nbody=converter)
g.particles.add_particles(p)
g.commit_particles()

E0 = g.kinetic_energy + g.potential_energy
print 'E0 =', E0

# Include mass loss.

mdotfac = 0.
for i in range(n):
    g.set_dmdt(i, -mdotfac*p[i].mass/(1000|units.Myr))

g.parameters.epsilon_squared = (0.01|units.parsec)**2
g.parameters.integrator = 5   			# 5th order symplectic
g.parameters.timestep = 0.005|units.Myr
g.parameters.timestep_parameter = 0.02	# timestep and timestep_parameter are
					# mutually exclusive -- no longer
                                        # symplectic if timestep_parameter > 0
print g.parameters

g.evolve_model(100.|units.Myr)

print 'time =', g.model_time.in_(units.Myr)
#print g.particles
E1 = g.kinetic_energy + g.potential_energy
print 'dE/E =', E1/E0 - 1

g.stop()
