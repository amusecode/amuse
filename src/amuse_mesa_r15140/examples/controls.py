import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse import datamodel

stellar_evolution = MESA(version='15140')

masses=[5.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

print(star.get_control('mesh_delta_coeff'))
star.set_control('mesh_delta_coeff',0.8)
print(star.get_control('mesh_delta_coeff'))

print(star.get_control('merge_if_dr_div_cs_too_small'))
star.set_control('merge_if_dr_div_cs_too_small',False)
print(star.get_control('merge_if_dr_div_cs_too_small'))

print(star.get_control('terminal_show_age_units'))
star.set_control('terminal_show_age_units','seconds')
print(star.get_control('terminal_show_age_units'))

print(star.get_control('solver_save_photo_call_number'))
star.set_control('solver_save_photo_call_number',3)
print(star.get_control('solver_save_photo_call_number'))


print(star.get_control('max_num_subcells'))
star.set_control('max_num_subcells',3)
print(star.get_control('max_num_subcells'))


# On their own changing an option in star_job does not do anything
# Once you've set all the star_job options you then want then call star_job_update()

print(star.get_star_job('min_x_for_keep'))
star.set_star_job('min_x_for_keep',10**-4)
print(star.get_star_job('min_x_for_keep'))

print(star.get_star_job('change_net'))
star.set_star_job('change_net',True)
print(star.get_star_job('change_net'))

print(star.get_star_job('new_net_name'))
star.set_star_job('new_net_name','approx21.net')
print(star.get_star_job('new_net_name'))

print(star.get_star_job('new_rates_preference'))
star.set_star_job('new_rates_preference',1)
print(star.get_star_job('new_rates_preference'))

star.star_job_update()


# EOS and KAP options

print(star.get_eos('mass_fraction_limit_for_PC'))
star.set_eos('mass_fraction_limit_for_PC',0.0)
print(star.get_eos('mass_fraction_limit_for_PC'))

print(star.get_kap('zbase'))
star.set_kap('zbase',0.001)
print(star.get_kap('zbase'))



print(star.get_history('species')) # Should be 21 now

# Change the mass of the tar
star.mass=10| units.MSun

print(star.mass)
print(star.get_history('star_mass'))



#Add extra energy to satr

print(star.get_extra_heat())
print(star.set_extra_heat(1.0,1))
print(star.get_mesa_value('extra_heat',1))

eng = np.zeros(star.get_number_of_zones())
eng[:] = 2.0
print(star.set_mesa_value_profile('extra_heat',eng))
print(star.get_mesa_value_profile('extra_heat'))

