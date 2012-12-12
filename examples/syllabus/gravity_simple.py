"""
   Simple routine for running a gravity code
"""
from amuse.lab import *

def main(N=10, W0=7.0, t_end=10, dt=1, filename="nbody.hdf5"):
    t_end = t_end | nbody_system.time
    dt = dt | t_end.unit
    bodies = new_king_model(N, W0)
    bodies.scale_to_standard()
    bodies.radius = 0 |  nbody_system.length
    gravity = ph4()
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)
    write_set_to_file(bodies.savepoint(0.0 | t_end.unit), filename, "hdf5")

    Etot_prev = Etot_init = gravity.kinetic_energy + gravity.potential_energy
    time = 0.0 | t_end.unit
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
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-t", dest="t_end", type="float", default = 1,
                      help="end time of the simulation [1] N-body units")
    result.add_option("-d", dest="dt", type="float", default = 0.1,
                      help="diagnostics time step [0.1] N-body units")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the potential (W0) [7.0]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

