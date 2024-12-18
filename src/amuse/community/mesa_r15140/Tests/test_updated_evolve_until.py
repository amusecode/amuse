from amuse.lab import *
import numpy as np
import matplotlib.pyplot as plt


# This script evolves a 10 Msun single star with MESA
# To test the updated evolve function (evolve_until in mesa_interface.f90)

stars = Particles(1)

stars[0].mass = 10 | units.MSun
stars[0].metallicity = 0.02


stellar_code = MESA(redirection="none")

evo_stars = stellar_code.particles.add_particles(stars)

# Standard MESA setup
stellar_code.particles[0].set_control('initial_z',0.02)
stellar_code.particles[0].set_kap('use_Type2_opacities',True)
stellar_code.particles[0].set_kap('Zbase', 0.02)


# To store time and timestep
times = np.array([])
time_step = np.array([])


time = 0 | units.Myr


# Subsequent calls to evolve function (the updated function is evolve_until in mesa_interface.f90)
# To check that the timestep can increase
# In the previous version it would be prevented from increasing when the evolve function was called
# multiple times in a row

# With the update, it increases freely up to the maximum allowed value, i.e. the value given as argument
# to evolve_for


while time < 0.02 | units.Myr:
	time = (stellar_code.particles[0].age)
	dt = stellar_code.particles[0].time_step

	times = np.append(times,time.value_in(units.kyr))
	time_step = np.append(time_step,dt.value_in(units.kyr))

	stellar_code.evolve_for(1,1|units.kyr)



while time < 0.06 | units.Myr:
	time = (stellar_code.particles[0].age)
	dt = stellar_code.particles[0].time_step

	times = np.append(times,time.value_in(units.kyr))
	time_step = np.append(time_step,dt.value_in(units.kyr))

	stellar_code.evolve_for(1,5|units.kyr)


while time < 0.12 | units.Myr:
	time = (stellar_code.particles[0].age)
	dt = stellar_code.particles[0].time_step

	times = np.append(times,time.value_in(units.kyr))
	time_step = np.append(time_step,dt.value_in(units.kyr))

	stellar_code.evolve_for(1,10|units.kyr)



# After the evolution, the timesteps are plotted as function of time. The three visible plateau
# correspond to the value given as argument to evolve_for (1 kyr, 5 kyrs, 10 kyrs),
# and indicate that the timestep was able to reach the maximum allowed value.

# When running the same script with the previous evolve_until subroutine, the time step is cut at
# the end of each call to evolve_until and can not increase freely to reach the maximum allowed value,
# which is why the plateau are not seen in this case.

plt.plot(times,time_step,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'd$t$ [kyr]',fontsize=13)
plt.show()
