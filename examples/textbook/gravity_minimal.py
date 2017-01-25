from amuse.lab import Hermite 
from amuse.lab import nbody_system
from amuse.lab import new_king_model

def main(N, W0, t_end):
    bodies = new_king_model(N, W0)
    bodies.scale_to_standard()

    gravity = Hermite()
    gravity.particles.add_particles(bodies)
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    gravity.evolve_model(t_end)

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    print "T=", gravity.get_time(), "M=", bodies.mass.sum(), 
    print "E= ", Etot, "Q= ", Ekin/Epot, "dE=", (Etot_init-Etot)/Etot

    gravity.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100, help="number of stars [%default]")
    result.add_option("-t", dest="t_end", unit=nbody_system.time, type="float", default = 1, help="end time of the simulation [%default] N-body units")
    result.add_option("-W", dest="W0", type="float", default = 7.0, help="Dimension-less depth of the potential (W0) [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
