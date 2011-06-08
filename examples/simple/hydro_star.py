from amuse.support.data.core import Particle
from amuse.support.units import units
from amuse.community.evtwin.interface import EVtwin
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.plot import sph_particles_plot, native_plot


def create_particles():
    star =  Particle()
    star.mass = 3.0 | units.MSun
    
    stellar_evolution = EVtwin()
    se_star = stellar_evolution.particles.add_particle(star)
    
    print "Evolving", star.mass, "star with", stellar_evolution.__class__.__name__, "up to", 100 | units.Myr
    stellar_evolution.evolve_model(100 | units.Myr)
    
    print "Creating SPH particles from the (1D) stellar evolution model"
    sph_particles = convert_stellar_model_to_SPH(
        se_star, 
        1000
    ).gas_particles
    stellar_evolution.stop()
    return sph_particles


if __name__ in ("__main__", "__plot__"):
    sph_particles = create_particles()
    native_plot.figure(figsize = (6, 6), dpi = 50)
    sph_particles_plot(sph_particles)
    native_plot.show()
