import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_next.interface import MESA
from amuse import datamodel

stellar_evolution = MESA()

# Set Z for all stars?
stellar_evolution.parameters.metallicity = 0.0142

masses=[1.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.pure_he_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step()

print(star.time_step)

star.time_step = 1. | units.yr

star.evolve_one_step()

print(star.time_step)

h1 = star.get_profile('h1')
print(h1[0],h1[-1])

he4 = star.get_profile('he4')
print(he4[0],he4[-1])
