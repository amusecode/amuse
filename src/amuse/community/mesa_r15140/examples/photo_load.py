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

print(star.time_step)

star.time_step = 1. | units.yr

star.evolve_one_step()

star.save_photo('photo')
star.save_model('model.mod')


star.evolve_one_step()

age = star.age

model_number = star.get_history('model_number')


masses=[1.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses,filename=['photo'])

stars = stellar_evolution.photo_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step() # Photos will give exactly the same answer compared to a run without a save/load

print(star.age, star.get_history('model_number'), star.mass)


masses=[1.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses,filename=['model.mod'])

stars = stellar_evolution.pre_built_stars.add_particles(stars)
star = stars[0]

star.evolve_one_step() # Models do not give exactly the same answer as a photo

print(star.age, star.get_history('model_number'), star.mass)