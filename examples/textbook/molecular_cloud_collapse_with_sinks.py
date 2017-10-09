from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver

import numpy
  
from matplotlib import pyplot 

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from hydrodynamics_class_SE import Hydro

def run_molecular_cloud(N=100, Mcloud=100. | units.MSun, Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
    gas=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube).result
    gas.name = "gas"

    rho_cloud = Mcloud/Rcloud**3
    tff = 0.5427/numpy.sqrt(constants.G*rho_cloud)
    print "Freefall timescale=", tff.in_(units.Myr)
    #omega = 1./tff
    #print "Omega=", omega.in_(units.yr**-1)
    #gas.add_spin(omega)

    stars = Particles(0)
    
    hydro = Hydro(Fi, gas, stars)

    rho_cloud = 3.*Mcloud/(4.*numpy.pi*Rcloud**3)
    print rho_cloud

#    dt = 0.1*tff
#    tend=40.0*tff
    dt = 0.02 | units.Myr
    tend= 10.0 | units.Myr
    dt_diag = 0.1 | units.Myr
    t_diag = 0 | units.Myr

    i=0
    E0 = 0.0
    time = 0.0 | units.Myr

    while time < tend:
        time += dt
        print "Evolve to time=", time.in_(units.Myr)
#        print "N=", len(sph.gas_particles), len(sph.dm_particles)

        Mtot = 0|units.MSun
        if len(hydro.star_particles)>0:
            print "Mass conservation: Slocal:", time.in_(units.Myr), hydro.gas_particles.mass.sum().in_(units.MSun), hydro.star_particles.mass.sum().in_(units.MSun), "sum=", (hydro.gas_particles.mass.sum()+hydro.star_particles.mass.sum()).in_(units.MSun)
            print "Mass conservation: Shydro:", time.in_(units.Myr), hydro.code.gas_particles.mass.sum().in_(units.MSun), hydro.code.dm_particles.mass.sum().in_(units.MSun), "sum=", (hydro.code.gas_particles.mass.sum()+hydro.code.dm_particles.mass.sum()).in_(units.MSun), "S=", hydro.stellar.particles.mass.sum().in_(units.MSun)
            Mtot = hydro.gas_particles.mass.sum()+hydro.star_particles.mass.sum()
        else:
            print "Mass conservation: local:", time.in_(units.Myr), hydro.gas_particles.mass.sum().in_(units.MSun) 
            print "Mass conservation: hydro:", time.in_(units.Myr), hydro.code.gas_particles.mass.sum().in_(units.MSun)
            Mtot = hydro.gas_particles.mass.sum()

        if Mtot<Mcloud-(1.e-5|units.MSun):
            print "Mass is not conserved:", Mtot.in_(units.MSun), Mcloud.in_(units.MSun)
            exit(-1)

        hydro.evolve_model(time)
        E = hydro.gas_particles.kinetic_energy()+hydro.gas_particles.potential_energy() + hydro.gas_particles.thermal_energy()
        E_th = hydro.gas_particles.thermal_energy()
        if i==0:
            E0 = E
        Eerr = (E-E0)/E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        print "maximal_density:",gas.rho.max().in_(units.MSun/units.parsec**3)

        hydro.print_diagnostics()
        if time>t_diag:
            t_diag += dt_diag
            hydro.write_set_to_file(i)
            i=i+1

    hydro.stop()
    return gas
  
if __name__ in ("__main__","__plot__"):

    numpy.random.seed(3141)
#    parts = run_molecular_cloud(10000, Mcloud=1000. | units.MSun, Rcloud=1.4 | units.parsec)
    parts = run_molecular_cloud(10000, Mcloud=1000. | units.MSun, Rcloud=3.0 | units.parsec)


    
