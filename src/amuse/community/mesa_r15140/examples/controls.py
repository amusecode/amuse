import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa_r15140.interface import MESA
from amuse import datamodel

stellar_evolution = MESA()

masses=[5.0] | units.MSun
stars = datamodel.Particles(len(masses), mass=masses)

stars = stellar_evolution.native_stars.add_particles(stars)
star = stars[0]

print(star.get_control_dble('mesh_delta_coeff'))
star.set_control_dble('mesh_delta_coeff',0.8)
print(star.get_control_dble('mesh_delta_coeff'))

print(star.get_control_logical('merge_if_dr_div_cs_too_small'))
star.set_control_logical('merge_if_dr_div_cs_too_small',False)
print(star.get_control_logical('merge_if_dr_div_cs_too_small'))

print(star.get_control_str('terminal_show_age_units'))
star.set_control_str('terminal_show_age_units','seconds')
print(star.get_control_str('terminal_show_age_units'))

print(star.get_control_int('max_num_subcells'))
star.set_control_int('max_num_subcells',3)
print(star.get_control_int('max_num_subcells'))


# On their own changing an option in star_job does not do anything
# Once you've set all the star_job options you then want then call star_job_update()

print(star.get_star_job_dble('min_x_for_keep'))
star.set_star_job_dble('min_x_for_keep',10**-4)
print(star.get_star_job_dble('min_x_for_keep'))

print(star.get_star_job_logical('change_net'))
star.set_star_job_logical('change_net',True)
print(star.get_star_job_logical('change_net'))

print(star.get_star_job_str('new_net_name'))
star.set_star_job_str('new_net_name','approx21.net')
print(star.get_star_job_str('new_net_name'))

print(star.get_star_job_int('new_rates_preference'))
star.set_star_job_int('new_rates_preference',1)
print(star.get_star_job_int('new_rates_preference'))


star.star_job_update()
print(star.get_history('species')) # Should be 21 now

# Change the mass of the tar
star.mass=10| units.MSun

print(star.mass)
print(star.get_history('star_mass'))