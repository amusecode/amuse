"""
   Minimalistic routine for running a hydrodynamics solver
"""
from amuse.lab import *
from amuse.io import store

def main(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun, t_end=1|units.Myr):
    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    gas = new_plummer_gas_model(N, convert_nbody=converter)
    gas.move_to_center()
    gas.vx = 100 | units.kms

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(gas)
    Etot_init = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    hydro.evolve_model(t_end)

    Ekin = hydro.kinetic_energy 
    Epot = hydro.potential_energy
    Eth = hydro.thermal_energy
    Etot = Ekin + Epot + Eth
    print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), 
    print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot
    print "pos-vel=", hydro.gas_particles.center_of_mass(), hydro.gas_particles.center_of_mass_velocity()

    hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 1|units.Myr,
                      help="end time of the simulation [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1|units.MSun,
                      help="Mass of molcular cloud [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="Rvir", type="float", default = 1|units.RSun,
                      help="Radius of cloud [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

