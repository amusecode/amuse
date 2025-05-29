import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse import datamodel

stellar_evolution = MESA(version='15140')

stellar_evolution.parameters.metallicity = 0.0142

masses=[1.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.pre_ms_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step()

star.evolve_for(4.6 | units.Gyr)

print(star.temperature,star.mass,star.radius,star.luminosity)

star.evolve_one_step()

print(star.temperature,star.get_history('Teff'),star.get_history('log_Teff'))

h1 = star.get_profile('h1')
print(h1[0],h1[-1])

# Note use the name from the history_columns.list file 'center h1' not the hsitory output 'center_h1'
print(star.get_history('center h1')) 
