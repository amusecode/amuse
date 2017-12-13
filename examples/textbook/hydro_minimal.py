###BOOKLISTSTART1###
from amuse.lab import *

def main(N, Mtot, Rvir, t_end):
    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    gas = new_plummer_gas_model(N, convert_nbody=converter)

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(gas)
    Etot_init = hydro.kinetic_energy \
              + hydro.potential_energy + hydro.thermal_energy
    hydro.evolve_model(t_end)
    write_set_to_file(hydro.particles, "hydro.h5", "hdf5")

    Ekin = hydro.kinetic_energy 
    Epot = hydro.potential_energy
    Eth = hydro.thermal_energy
    Etot = Ekin + Epot + Eth
    Q = (Ekin+Eth)/Epot
    dE = (Etot_init-Etot)/Etot
    com = hydro.gas_particles.center_of_mass()
    print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), 
    print "E= ", Etot, "Q= ", Q, "dE=", dE, "CoM=", com.in_(units.RSun)

    hydro.stop()
###BOOKLISTSTOP1###
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of gas particles [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 6|units.hour,
                      help="end time of the simulation [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1|units.MSun,
                      help="Mass of the cloud [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="Rvir", type="float", default = 1|units.RSun,
                      help="Radius of the cloud [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

