import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_next.interface import MESA
from amuse import datamodel

stellar_evolution = MESA()

masses=[5.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)

# Show how to make a star based on an importted model
# For convenience we just use another star made with MESA
# but this could be from another code

# Lets evolve a star first
stars[0].evolve_for(0.05 | units.Gyr)


# Make new one based on what we just created
new_star = stellar_evolution.new_particle_from_model(stars[0], 0.0|units.Myr)


# Check they look similiar, they wont be identical as mesa attemtps to find a matching solution
# within some tolerances (even when importing an exisiting mesa model), this will be worse
# the more the EOSs and other microphysics differs between models. 

print(new_star.mass,stars[0].mass)

print(new_star.temperature,stars[0].temperature)

print(new_star.get_profile('h1')[0],stars[0].get_profile('h1')[0])

#print(new_star.age,stars[0].age)


new_star.evolve_one_step()
#print(new_star.age)
