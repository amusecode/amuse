import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse import datamodel

stellar_evolution = MESA(version='15140')

masses=[2.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step()

print(star.get_number_of_species())
print(star.get_name_of_species(1))

print(star.get_masses_of_species())

print(star.get_mass_fraction_of_species_at_zone(1,1))

print(star.get_id_of_species('h1'))
