"""
Generates an isochrone of a cluster of stars in a Hertzsprung-Russell diagram
"""
from matplotlib import pyplot
from amuse.units import units
from amuse.datamodel import Particles
from amuse.ic.brokenimf import new_scalo_mass_distribution
from amuse.ext.particles_with_color import new_particles_with_blackbody_color
from amuse.community.seba.interface import SeBa

def simulate_stellar_evolution(particles, endtime):
    stellar_evolution = SeBa()
    stellar_evolution.particles.add_particles(particles)
    from_code_to_model = stellar_evolution.particles.new_channel_to(particles)
    print "Evolving {0} stars using {1} up to {2}.".format(len(particles), stellar_evolution.__class__.__name__, endtime)
    stellar_evolution.evolve_model(endtime)
    from_code_to_model.copy_all_attributes()
    stellar_evolution.stop()

def plot_isochrone(particles):
    particles = new_particles_with_blackbody_color(particles)
    pyplot.figure(figsize = (7, 8))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    
    pyplot.scatter(
        particles.temperature.value_in(units.K), 
        particles.luminosity.value_in(units.LSun), 
        s=particles.radius.maximum(0.1|units.RSun).value_in(units.RSun)*100, 
        c=particles.color, 
        edgecolors = "none", 
    )
    pyplot.xlabel('Effective Temperature (K)')
    pyplot.ylabel('Luminosity (L$_{\odot}$)')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlim(20e3, 2e3)
    pyplot.ylim(1.0e-3, 1.0e3)
    pyplot.gca().set_axis_bgcolor('#808080')
    pyplot.show()

if __name__ in ('__main__', '__plot__'):
    particles = Particles(mass=new_scalo_mass_distribution(4000))
    simulate_stellar_evolution(particles, 0.5 | units.Gyr)
    plot_isochrone(particles)
