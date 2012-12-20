"""
Convert a 1D Henhey stellar model into a spherical particle distribution
"""
from amuse.units import units
from amuse.community.evtwin.interface import EVtwin
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.plot import sph_particles_plot, native_plot
from amuse.datamodel import Particle

def convert_star_to_hydro_model(M, t_end):

    star =  Particle(mass=M)
    stellar_evolution = EVtwin()
    se_star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.evolve_model(t_end)
    
    sph_particles = convert_stellar_model_to_SPH(se_star, 10000).gas_particles
    stellar_evolution.stop()
    return sph_particles

if __name__ in ("__main__", "__plot__"):
    sph_particles = convert_star_to_hydro_model(2.0|units.MSun, 1100|units.Myr)
    native_plot.figure(figsize = (10, 10), dpi = 50)
    sph_particles_plot(sph_particles)
    native_plot.show()
