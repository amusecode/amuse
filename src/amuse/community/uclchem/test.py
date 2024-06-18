from interface import Uclchem
from amuse.datamodel import Particles
from amuse.units import units

chem = Uclchem()
particles = Particles(1)
particles.number_density = 1.0 | units.cm**-3
particles.temperature = 20 | units.K
particles.ionrate = 10**-16 | units.s**-1

print(particles)

dt = 1e5|units.yr
chem.particles.add_particle(particles)
channel = chem.particles.new_channel_to(particles)
t=0|units.yr
index_H = chem.get_index_of_species('H')
for i in range(10):
    t += dt
    chem.evolve_model(t)
    channel.copy()
    print(chem.particles.abundances[0][index_H])
print(particles)
    
