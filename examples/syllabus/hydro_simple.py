"""
   Simple routine for running a hydrodynamics solver
"""
from amuse.lab import *
from amuse.io import store

def main(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun, 
         t_end=1|units.Myr, n_steps=10):

    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(bodies)
    Etot_init = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)
    write_set_to_file(bodies.savepoint(0.0 | t_end.unit), "hydro.hdf5", "hdf5")

    time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)
    while time < t_end:
        time += dt
        hydro.evolve_model(time)
        hydro_to_framework.copy()
        write_set_to_file(bodies.savepoint(time), "hydro.hdf5", "hdf5")

        Ekin = hydro.kinetic_energy 
        Epot = hydro.potential_energy
        Eth = hydro.thermal_energy
        Etot = Ekin + Epot + Eth
        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), 
        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot
    hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-n", dest="n_steps", type="int",default = 10,
                      help="number of steps [10]")
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

