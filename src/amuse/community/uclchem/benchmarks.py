from amuse.datamodel import Particles
from amuse.units import units, nbody_system, constants
from amuse.community.uclchem.interface import UCLchem
from amuse.community.fi.interface import Fi
from amuse.ext.molecular_cloud import molecular_cloud

import matplotlib.pyplot as plt
import numpy as np

def run_benchmark(model='static', n_init=1e4|units.cm**-3, n_end=1e6|units.cm**-3, filename='static_cloud_benchmark.dat'):
    R_cloud = 0.05|units.pc
    M_cloud = 4/3*constants.pi*R_cloud**3*units.amu*n_init
    T_cloud = 10.0|units.K
    cr_ion_cloud = 1|units.cr_ion

    conv = nbody_system.nbody_to_si(R_cloud, M_cloud)
#cloud = molecular_cloud(targetN=1000, convert_nbody=conv).result

    cloud = Particles(1)
    cloud.number_density = 1.001*n_init
    cloud.temperature = T_cloud
    cloud.ionrate = cr_ion_cloud

    chem = UCLchem()
    chem.out_species = ["OH", "OCS", "CO", "CS", "CH3OH"]
    chem.particles.add_particles(cloud)
    chem_channel = chem.particles.new_channel_to(cloud)
    channel_to_chem = cloud.new_channel_to(chem.particles)

    t_end = 5.0e6|units.yr
    t = 0|units.yr
    #times = np.logspace(-7,np.log10(5e6),200)
    first, last = chem.get_firstlast_abundance()
    species = 'time,density'
    for i in range(last):
        species = species + ', ' + chem.get_name_of_species(i)
    print(chem.particles.number_density)
    results = np.insert(chem.particles.abundances[0], 0,[t.value_in(units.yr),chem.particles.number_density[0].value_in(units.cm**-3)])
    print(results)
    if model =='static':
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
        results = np.vstack((results, np.insert(chem.particles.abundances[0], 0, [t.value_in(units.yr),chem.particles.number_density[0].value_in(units.cm**-3)])))
        
    elif model =='freefall':
        t_old = t
        while chem.particles.number_density < n_end:
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
            cloud.number_density = freefall_eq(chem.particles.number_density[0], n_init, t-t_old)
            channel_to_chem.copy()
            print(chem.particles.number_density)
            chem.evolve_model(t)
            chem_channel.copy()
            results = np.vstack((results, np.insert(chem.particles.abundances[0], 0, [t.value_in(units.yr),chem.particles.number_density[0].value_in(units.cm**-3)])))
            t_old = t
            

    np.savetxt(filename,results,header=species, delimiter = ',', comments=' ')
    chem.stop()

def freefall_eq(density, initial_density, dt):
    dndt = (density**4/initial_density)**(1/3)*(24*np.pi*constants.G*constants.u*initial_density*((density/initial_density)**(1/3)-1))**(1/2)
    dn = dndt*dt
    return density+dn



run_benchmark(model='freefall', n_init=1e2|units.cm**-3,filename='freefall_benchmark.dat')