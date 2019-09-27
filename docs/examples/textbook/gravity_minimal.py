import time
from amuse.lab import Hermite 
from amuse.lab import nbody_system
from amuse.lab import new_king_model

def gravity_minimal(N, W0, t_end):
    bodies = new_king_model(N, W0)
    bodies.scale_to_standard()

    gravity = Hermite()
    gravity.particles.add_particles(bodies)
    Etot_init = gravity.kinetic_energy + gravity.potential_energy

    start_time = time.time()
    gravity.evolve_model(t_end)
    dtime = time.time() - start_time 

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy

    Etot = Ekin + Epot
    dE = (Etot_init-Etot)/Etot
    print "T =", gravity.get_time(), " CPU time:", dtime, "[s]"
    print "M =", bodies.mass.sum(), " E = ", Etot, " Q = ", -Ekin/Epot
    print "dE =", dE
    
    gravity.stop()
    
if __name__ in ('__main__'):
    N = 100
    W0 = 7.0
    t_end = 1 | nbody_system.time
    gravity_minimal(N, W0, t_end)

