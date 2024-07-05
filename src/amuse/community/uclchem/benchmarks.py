from amuse.datamodel import Particles
from amuse.units import units, nbody_system, constants
from amuse.community.uclchem.interface import UCLchem
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

chem = UCLchem()
chem.out_species = ["OH", "OCS", "CO", "CS", "CH3OH"]
chem.particles.add_particles(static_cloud)
chem_channel = chem.particles.new_channel_to(static_cloud)

t_end = 5.0e6|units.yr
t = 0|units.yr
times = np.logspace(-7,np.log10(5e6),200)
first, last = chem.get_firstlast_abundance()
species = 'time'
for i in range(last):
    species = species + ', ' + chem.get_name_of_species(i)

results = np.insert(chem.particles.abundances[0], 0,t.value_in(units.yr))

while t.value_in(units.yr) < t_end.value_in(units.yr):
    if t >= 1.0e6|units.yr:
        t += 1.0e5|units.yr
    elif t >= 1.0e5|units.yr:
        t += 1.0e4|units.yr
    elif t >= 1.0e4|units.yr:
        t += 1.0e3|units.yr
    elif t >= 1.0e2|units.yr:
        t += 100|units.yr
    elif t > 0|units.yr:
        t *= 10
    else:
        t = 1.0e-7|units.yr
# for t in times:
#     print(t)
#     t = t|units.yr
    chem.evolve_model(t)
    chem_channel.copy()
    results = np.vstack((results, np.insert(chem.particles.abundances[0], 0, t.value_in(units.yr))))

np.savetxt('static_cloud_benchmark.dat',results,header=species, delimiter = ',', comments=' ')
chem.stop()