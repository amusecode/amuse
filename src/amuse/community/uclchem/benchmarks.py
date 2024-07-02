from amuse.datamodel import Particles
from amuse.units import units, nbody_system, constants
from amuse.community.uclchem.interface import Uclchem
from amuse.community.fi.interface import Fi
from amuse.ext.molecular_cloud import molecular_cloud

import matplotlib.pyplot as plt
import numpy as np

R_cloud = 0.05|units.pc
n_cloud = 1e4|units.cm**-3
M_cloud = 4/3*constants.pi*R_cloud**3*units.amu*n_cloud
T_cloud = 10.0|units.K
cr_ion_cloud = 1|units.cr_ion

conv = nbody_system.nbody_to_si(R_cloud, M_cloud)
#cloud = molecular_cloud(targetN=1000, convert_nbody=conv).result

static_cloud = Particles(1)
static_cloud.number_density = n_cloud
static_cloud.temperature = T_cloud
static_cloud.ionrate = cr_ion_cloud

chem = Uclchem()
chem.out_species = ["OH", "OCS", "CO", "CS", "CH3OH"]
chem.particles.add_particles(static_cloud)
chem_channel = chem.particles.new_channel_to(static_cloud)

t_end = 5.0e6|units.yr
#dt = np.logspace(0,np.log10(t_end.value_in(units.yr)), num = 200)
#dt = np.linspace(0,5e6, num=200)
t = 0|units.yr
times = []
index_H = chem.get_index_of_species('H')
index_H2 = chem.get_index_of_species('H2')
index_iceH = chem.get_index_of_species('#H')
index_bulkH = chem.get_index_of_species('@H')
index_H2O = chem.get_index_of_species('H2O')
index_iceH2O = chem.get_index_of_species('#H2O')
index_bulkH2O = chem.get_index_of_species('@H2O')
index_CO = chem.get_index_of_species('CO')
index_iceCO = chem.get_index_of_species('#CO')
index_bulkCO = chem.get_index_of_species('@CO')
index_iceCH3OH = chem.get_index_of_species('#CH3OH')
index_bulkCH3OH = chem.get_index_of_species('@CH3OH')
index_CH3OH = chem.get_index_of_species('CH3OH')
indices = np.array([index_H,index_H2,index_iceH,index_bulkH, index_H2O,index_iceH2O,index_bulkH2O,index_CO,index_bulkCO,index_iceCH3OH,index_bulkCH3OH,index_CH3OH])
abundances = np.zeros(len(indices))
names = ['H','H2','$H','H2O','$H2O','CO','$CO','$CH3OH','CH3OH']
print(indices)


while t.value_in(units.yr) < t_end.value_in(units.yr):
    if t > 1.0e6|units.yr:
        t += 1.0e5|units.yr
    elif t > 1.0e5|units.yr:
        t += 1.0e4|units.yr
    elif t > 1.0e4|units.yr:
        t += 1.0e3|units.yr
    elif t > 1.0e3|units.yr:
        t += 100|units.yr
    elif t > 0|units.yr:
        t *= 10
    else:
        t = 1.0e-7|units.yr
    chem.evolve_model(t)
    abundances = np.vstack((abundances,chem.particles.abundances[0][indices]))
    chem_channel.copy()
    times.append(t.value_in(units.yr))

# for i in range(len(indices)):
#     if i in [0,1,3,5,8]:
#         plt.plot(dt,abundances[1:,i],label=names[i])
#     else: 
#         plt.plot(dt,abundances[1:,i], linestyle='dashed',label=names[i])
print(chem.particles.abundances[0][indices])
plt.plot(times,abundances[1:,0], label='H')
plt.plot(times,abundances[1:,1], label='H2')
plt.plot(times,abundances[1:,2]+abundances[1:,3], linestyle='dashed',label='$H')
plt.plot(times,abundances[1:,4], label='H2O')
plt.plot(times,abundances[1:,5]+abundances[1:,6],linestyle='dashed', label='$H2O')
plt.plot(times, abundances[1:,7],label='CO')
plt.plot(times,abundances[1:,8]+abundances[1:,9],linestyle='dashed',label='$CO')
plt.plot(times,abundances[1:,10]+abundances[1:,11],linestyle='dashed',label='$CH3OH')
plt.plot(times,abundances[1:,11], label='CH3OH')
plt.title('Static Cloud')
plt.xlabel('Time (yr)') 
plt.ylabel('Abundances')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('static_cloud.pdf')
plt.show()
plt.close()