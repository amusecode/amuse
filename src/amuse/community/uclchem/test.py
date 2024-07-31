from interface import UCLchem
from amuse.datamodel import Particles
from amuse.units import units

chem = UCLchem()
particles = Particles(2)
particles[0].number_density = 1.0e2 | units.cm**-3
particles[0].temperature = 20 | units.K
particles[0].ionrate = 10**-16 | units.s**-1
particles[0].radfield = 1.5 | units.habing

particles[1].number_density = 1.0e3 | units.cm**-3
particles[1].temperature = 10 | units.K
particles[1].ionrate = 5*10**-16 | units.s**-1
particles[1].radfield = 3 | units.habing
print(particles)
chem.out_species = ['H','H2', '@H2']

print('species', chem.out_species)
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
    
