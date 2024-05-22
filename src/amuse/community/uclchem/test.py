from interface import Uclchem
from amuse.datamodel import Particles
from amuse.units import units
species='H H2'
dict = "{'outspecies': 2}"
chem = Uclchem()
particles = Particles(1)
particles.dens = 1.0 | units.g * units.cm**-3
particles.temperature = 20 | units.K
particles.ionrate = 10**-16 | units.s**-1
print(particles)
chem.particles.add_particles(particles)
chem.particle_test()
#print(chem.sim_cloud(outSpeciesIn=species, dictionary=dict))
#print(test)