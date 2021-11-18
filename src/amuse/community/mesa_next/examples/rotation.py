import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_next.interface import MESA
from amuse import datamodel

stellar_evolution = MESA()

masses=[1.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

# Turn on rotation
star.set_star_job('change_rotation_flag',True)
star.set_star_job('new_rotation_flag',True)

# Many options for initialising the rotation rate
star.set_star_job('set_surface_rotation_v',True)
star.set_star_job('new_surface_rotation_v',100.0) # sets surface velocity in km/sec

# Now change the model to turn on rotation
star.star_job_update()

# Turn options back off
star.set_star_job_logical('change_rotation_flag',False)


# Evovle star until?
star.evolve_for(1.0 | units.Gyr)


irot = star.get_profile('irot') # ! specific moment of inertia at cell boundary
mass = star.get_profile('mass') # m/Msun. mass coordinate of outer boundary of cell.

