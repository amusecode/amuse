import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from amuse.io import read_set_from_file
from amuse.units import units
from amuse.community.uclchem.interface import UCLchem
from plotting import HydroPlotter

uclchem_data = pd.read_table('phase1-full.dat', skiprows=2, sep=',')
uclchem_data = uclchem_data.rename(columns=lambda x: x.strip())

static_cloud_data =pd.read_table('freefall_benchmark_2.dat', sep=',')
static_cloud_data = static_cloud_data.rename(columns=lambda x: x.strip())

plt.plot(static_cloud_data['time'],static_cloud_data['H'], color='b', label='H')
plt.plot(static_cloud_data['time'],static_cloud_data['H2'], color='r',label='H2')
plt.plot(static_cloud_data['time'],static_cloud_data['@H']+static_cloud_data['#H'],color='g',label='$H')
plt.plot(static_cloud_data['time'],static_cloud_data['H2O'],color='m', label='H2O')
plt.plot(static_cloud_data['time'],static_cloud_data['@H2O']+static_cloud_data['#H2O'],color='k', label='$H2O')
plt.plot(static_cloud_data['time'],static_cloud_data['CO'],color='c',label='CO')
plt.plot(static_cloud_data['time'],static_cloud_data['@CO']+static_cloud_data['#CO'],color='purple',label='$CO')
plt.plot(static_cloud_data['time'],static_cloud_data['@CH3OH']+static_cloud_data['#CH3OH'],color='grey',label='$CH3OH')
plt.plot(static_cloud_data['time'],static_cloud_data['CH3OH'],color='brown', label='CH3OH')

plt.plot(uclchem_data['Time'], uclchem_data['H'],linestyle='dashed',color='b')
plt.plot(uclchem_data['Time'], uclchem_data['H2'],linestyle='dashed',color='r')
plt.plot(uclchem_data['Time'], uclchem_data['@H']+uclchem_data['#H'],color='g',linestyle='dashed')
plt.plot(uclchem_data['Time'], uclchem_data['H2O'],linestyle='dashed',color='m')
plt.plot(uclchem_data['Time'], uclchem_data['@H2O']+uclchem_data['#H2O'],color='k',linestyle='dashed')
plt.plot(uclchem_data['Time'], uclchem_data['CO'],linestyle='dashed', color='c')
plt.plot(uclchem_data['Time'], uclchem_data['#CO']+uclchem_data['@CO'],color='purple',linestyle='dashed')
plt.plot(uclchem_data['Time'], uclchem_data['@CH3OH']+uclchem_data['#CH3OH'],color='grey',linestyle='dashed')
plt.plot(uclchem_data['Time'], uclchem_data['CH3OH'],linestyle='dashed', color='brown')

plt.title('Cloud in Freefall')
plt.xlabel('Time (yr)') 
plt.ylabel('Abundances')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlim(10,5e6)
plt.ylim(1e-15)
plt.savefig('freefall_zoom_2.pdf')
plt.show()
plt.close()

plt.plot(static_cloud_data['time'], static_cloud_data['density'])
plt.plot(uclchem_data['Time'], uclchem_data['Density'],linestyle='dashed')
plt.show()

CO_index = UCLchem().get_index_of_species('CO')
dt = 0.005
res = 500
plotter = HydroPlotter(res, use_torch=False)
for i in range(1,26):
    particles = read_set_from_file("mol_cloud_uclchem_{}.txt".format(i),format='amuse')
    plotter.add_gas_particles(particles)
    plt.plot(particles.density.value_in(units.g*units.cm**-3),particles.abundances[:,CO_index], '.',label='CO')
    plt.xlabel('Density (g/cm^3)')
    plt.ylabel('abundance')
    plt.title('abundances T={}Myr'.format(i*dt))
    plt.legend()
    plt.savefig('density-abundance_{}.pdf'.format(i))
    plt.show()

    plotter.plot_projection(2.5,params_to_plot=['density','abundance'],species_index=CO_index,use_gaussian_kernel=True, save_fig="projection_{}.pdf".format(i))
    plotter.clear_all_particles()