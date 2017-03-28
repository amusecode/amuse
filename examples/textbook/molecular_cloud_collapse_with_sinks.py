from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver

import numpy
  
from matplotlib import pyplot 

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from hydrodynamics_class import Hydro

def run_molecular_cloud(N=100, Mcloud=100. | units.MSun, Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
    gas=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube, seed=100).result
    gas.name = "gas"
    
    hydro = Hydro(Fi, gas)

    rho_cloud = 3.*Mcloud/(4.*numpy.pi*Rcloud**3)
    print rho_cloud
    tff = 0.5427/numpy.sqrt(constants.G*rho_cloud)
    print "Freefall timescale=", tff.in_(units.Myr)

    dt = 0.1*tff
    tend=4.0*tff

    i=0
    L=6
    E0 = 0.0
    time = 0.0 | units.Myr

    while time < tend:
        time += dt
        print "Evolve to time=", time.in_(units.Myr)
#        print "N=", len(sph.gas_particles), len(sph.dm_particles)
#        print "Masses of dm particles:", sph.dm_particles.mass.in_(units.MSun)

        hydro.evolve_model(time)
        E = hydro.gas_particles.kinetic_energy()+hydro.gas_particles.potential_energy() + hydro.gas_particles.thermal_energy()
        E_th = hydro.gas_particles.thermal_energy()
        if i==0:
            E0 = E
        Eerr = (E-E0)/E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        print "maximal_density:",gas.rho.max().in_(units.MSun/units.parsec**3)

        hydro.print_diagnostics()
        hydro.write_set_to_file(i)

        i=i+1

    hydro.stop()
    return gas
  
if __name__ in ("__main__","__plot__"):

    parts = run_molecular_cloud(1000, Mcloud=100. | units.MSun, Rcloud=3.0 | units.parsec)


    
