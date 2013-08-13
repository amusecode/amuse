from amuse.lab import *
from interface import nbody6xx
from amuse.units import nbody_system
inst = nbody6xx(redirection="none")

inst.initialize_code()
inst.particles.mass = 1 | nbody_system.mass

inst.evolve_model(2|nbody_system.time)

inst.get_velocity(1)
inst.get_position(1)
inst.get_number_of_particles()
inst.get_total_mass()
inst.get_total_radius()

print inst.particles
inst.particles.add_particles(new_plummer_model(10))
p=new_plummer_model(10)
inst.new_particle(p[0].mass,p[0].x,p[0].y,p[0].z,p[0].vx,p[0].vy,p[0].vz,p[0].radius)
