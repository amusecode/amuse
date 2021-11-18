import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_next.interface import MESA
from amuse import datamodel

stellar_evolution = MESA()

masses=[2.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step()

print(star.time_step)

star.time_step = 1. | units.yr

star.evolve_one_step()

print(star.get_history('time_step')) # Previous step

print(star.time_step) # Next step

