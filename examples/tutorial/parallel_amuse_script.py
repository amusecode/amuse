import time
import numpy
from amuse.lab import Huayno
from amuse.lab import Hermite
from amuse.lab import nbody_system
from amuse.lab import new_king_model
from matplotlib import pyplot

def gravity_minimal(bodies, t_end, nproc):

    gravity = Hermite(number_of_workers=nproc)
    gravity.particles.add_particles(bodies)
    Etot_init = gravity.kinetic_energy + gravity.potential_energy

    start_time = time.time()
    gravity.evolve_model(t_end)
    dtime = time.time() - start_time 

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy

    Etot = Ekin + Epot
    dE = (Etot_init-Etot)/Etot

    print()
    print("T =", gravity.get_time(), " CPU time:", dtime, "[s]")
    print("M =", bodies.mass.sum(), " E = ", Etot, " Q = ", -Ekin/Epot)
    print("dE =", dE)
    
    gravity.stop()
    return dtime
    
if __name__ in ('__main__'):
    N = 1024
    W0 = 7.0
    t_end = 0.1 | nbody_system.time
    bodies = new_king_model(N, W0)
    bodies.scale_to_standard()

    nproc= 6
    proc = numpy.arange(1, nproc+1, 1)
    tcpu = []
    for npi in proc:
        tcpu.append(gravity_minimal(bodies, t_end, npi))

    pyplot.scatter(proc, tcpu)
    pyplot.xlabel("n proc")
    pyplot.ylabel("CPU time [s]")
    pyplot.savefig("fig_parallel_performance_N1k_Hermite.pdf")

