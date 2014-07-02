from amuse.lab import *

def main(Ncl, rcl, W0, t_end, n_steps):

    masses = new_salpeter_mass_distribution(Ncl, 0.1|units.MSun, 100|units.MSun)
    converter = nbody_system.nbody_to_si(masses.sum(), rcl)
    bodies = new_king_model(Ncl, W0, convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(convert_nbody=converter)

    gravity = BHTree(converter)
    gravity.particles.add_particles(bodies)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)
    
    filename="nbody.hdf5"
    write_set_to_file(bodies.savepoint(0.0 | t_end.unit), filename, "hdf5", append_to_file=False)

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = zero
    dt = t_end/float(n_steps)
    while time < t_end:
        time += dt

        gravity.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        write_set_to_file(bodies.savepoint(time), filename, "hdf5")


        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, "M=", bodies.mass.sum(), 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

    gravity.stop()
    
def new_option_parser():
    
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Ncl", type="int",default = 100,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 1|units.Gyr,
                      help="end time of the simulation [%default]")
    result.add_option("-n", dest="n_steps", type="float", default = 100,
                      help="number of diagnostics output steps [%default]")
    result.add_option("-r", unit=units.parsec,
                      dest="rcl", type="float", default = 10|units.parsec,
                      help="cluster half-mass radius [%default]")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.RSun, units.yr], 
                          precision = 4, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

