from amuse.lab import Hermite 
from amuse.lab import nbody_system
from amuse.lab import new_plummer_model

def gravity_minimal(N, t_end):
    stars = new_plummer_model(N)

    gravity = Hermite()
    gravity.particles.add_particles(stars)

    initial_total_energy = gravity.kinetic_energy + gravity.potential_energy

    gravity.evolve_model(t_end)

    total_energy = gravity.kinetic_energy + gravity.potential_energy
    
    energy_error = (initial_total_energy-total_energy)/initial_total_energy
    
    print("model time =", gravity.model_time)
    print("total mass = ", stars.total_mass(), "total energy = ", total_energy)
    print("relative energy error =", energy_error)
        
if __name__ in ('__main__'):
    N = 100
    t_end = 1 | nbody_system.time
    gravity_minimal(N, t_end)

