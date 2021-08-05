import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_r15140.interface import MESA
from amuse import datamodel

stellar_evolution = MESA(gyre_in='./gyre.in')

masses=[1.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

star.time_step = 100. | units.yr

star.evolve_one_step()

print(star.get_history('time_step')) # Previous step

star.get_gyre()