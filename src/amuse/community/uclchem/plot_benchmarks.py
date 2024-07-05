import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

uclchem_data = pd.read_table('static-full.dat', skiprows=2, sep=',')
uclchem_data = uclchem_data.rename(columns=lambda x: x.strip())

static_cloud_data =pd.read_table('static_cloud_benchmark.dat', sep=',')
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

plt.title('Static Cloud')
plt.xlabel('Time (yr)') 
plt.ylabel('Abundances')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlim(10,5e6)
plt.ylim(1e-15)
plt.savefig('static_cloud_zoom.pdf')
plt.show()
plt.close()