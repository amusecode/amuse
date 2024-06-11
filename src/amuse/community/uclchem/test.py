from interface import Uclchem
from amuse.datamodel import Particles
from amuse.units import units
species='H H2'
dict = "{'outspecies': 2}"
chem = Uclchem()
particles = Particles(1)
particles.number_density = 1.0 | units.cm**-3
particles.temperature = 20 | units.K
particles.ionrate = 10**-16 | units.s**-1
chem.particles.add_particle(particles)
print(chem.particles.temperature)
chem.evolve_model(1e6|units.yr)
print(chem.particles.abundances)
