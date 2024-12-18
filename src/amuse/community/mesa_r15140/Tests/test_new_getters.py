from amuse.lab import *
import numpy as np
import matplotlib.pyplot as plt


# This script evolves a 10 Msun single star with MESA
# To test the getters added to the mesa interface

stars = Particles(1) 

stars[0].mass = 10 | units.MSun
stars[0].metallicity = 0.02


stellar_code = MESA(redirection="none")

evo_stars = stellar_code.particles.add_particles(stars)

# Standard MESA setup
stellar_code.particles[0].set_control('initial_z',0.02)
stellar_code.particles[0].set_kap('use_Type2_opacities',True)
stellar_code.particles[0].set_kap('Zbase', 0.02)

# Add wind mass loss
stellar_code.particles[0].set_control('hot_wind_scheme','Dutch')
stellar_code.particles[0].set_control('Dutch_scaling_factor',1.0)

# To store data
times = np.array([])

core_radius = np.array([])
conv_envelope_mass = np.array([])
conv_envelope_radius = np.array([])
wind_mass_loss_rate = np.array([])
apsidal_motion_constant = np.array([])
gyration_radius = np.array([])




time = 0 | units.Myr


# Evolve and store the data retrieved from the getters
while time < 0.1 | units.Myr:
	stellar_code.evolve_for(1,20|units.kyr)
	time = (stellar_code.particles[0].age)
	
	core_radius_temp = stellar_code.particles[0].core_radius
	conv_envelope_mass_temp = stellar_code.particles[0].convective_envelope_mass
	conv_envelope_radius_temp = stellar_code.particles[0].convective_envelope_radius
	wind_mass_loss_rate_temp = stellar_code.particles[0].wind_mass_loss_rate
	apsidal_motion_constant_temp = stellar_code.particles[0].apsidal_motion_constant
	gyration_radius_temp = stellar_code.particles[0].gyration_radius


	times = np.append(times,time.value_in(units.kyr))
	
	core_radius = np.append(core_radius,core_radius_temp.value_in(units.RSun))
	conv_envelope_mass = np.append(conv_envelope_mass,conv_envelope_mass_temp.value_in(units.MSun))
	conv_envelope_radius = np.append(conv_envelope_radius,conv_envelope_radius_temp.value_in(units.RSun))
	wind_mass_loss_rate = np.append(wind_mass_loss_rate,wind_mass_loss_rate_temp.value_in(units.MSun / units.yr))
	apsidal_motion_constant = np.append(apsidal_motion_constant,apsidal_motion_constant_temp)
	gyration_radius = np.append(gyration_radius,gyration_radius_temp)



# After the evolution, plot the data as function of time

plt.plot(times,core_radius,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'$R_{core}$ [R$_\odot$]',fontsize=13)
plt.show()

plt.plot(times,conv_envelope_mass,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'$M_{env}$ [M$_\odot$]',fontsize=13)
plt.show()

plt.plot(times,conv_envelope_radius,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'$R_{env}$ [R$_\odot$]',fontsize=13)
plt.show()

plt.plot(times,wind_mass_loss_rate,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'$\dot M$ [M$_\odot$/yr]',fontsize=13)
plt.show()

plt.plot(times,apsidal_motion_constant,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'$k_2$',fontsize=13)
plt.show()

plt.plot(times,gyration_radius,linewidth=2)
plt.xlabel(r'$t$ [kyr]',fontsize=13)
plt.ylabel(r'$r_g$',fontsize=13)
plt.show()