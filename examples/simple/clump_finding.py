import numpy 
"""
Generates a random distribution of particles and uses Hop to determine clumps
in the particle set.
"""
from numpy import random
from matplotlib import pyplot

from amuse.units import units
from amuse.units import nbody_system
from amuse.community.hop.interface import Hop

from amuse.datamodel.particles import Particles
from amuse.ic.salpeter import new_salpeter_mass_distribution

def new_cluster(number_of_stars = 1000):
    masses = new_salpeter_mass_distribution(
        number_of_stars, 
        mass_min = 0.1 | units.MSun,
        mass_max = 125.0 | units.MSun, 
        alpha = -2.35
    )
    
    particles = Particles(number_of_stars)
    particles.mass = masses
    particles.x = units.parsec(random.gamma(2.0, 1.0, number_of_stars))
    particles.y = units.parsec(random.gamma(1.0, 1.0, number_of_stars))
    particles.z = units.parsec(random.random(number_of_stars))
    
    return particles


def find_clumps(particles, unit_converter):
    
    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()
    hop.do_hop()
    
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    
    hop.stop()
    
    return result
    
def plot_clumps(groups, total_mass):
    number_of_particles_in_group = []
    fraction_of_mass_in_group =  []

    for group in groups:
        number_of_particles_in_group.append(len(group))
        fraction = (group.mass.sum()/total_mass)
        fraction_of_mass_in_group.append(fraction)
    
    figure = pyplot.figure(figsize= (12,6))
    
    
    subplot = figure.add_subplot(1, 2, 1)
    
    colormap = pyplot.cm.Paired
    for index, group in enumerate(groups):
        color = colormap(1.0 * index / len(groups))
        subplot.scatter(
            group.x.value_in(units.parsec),
            group.y.value_in(units.parsec),
            s = group.mass.value_in(units.MSun),
            edgecolors = color,
            facecolors = color
        )
    
    subplot.set_xlim(0,1)
    subplot.set_ylim(0,1)
    subplot.set_xlabel('x (parsec)')
    subplot.set_ylabel('y (parsec)')
    
    subplot = figure.add_subplot(1, 2, 2)
        
    subplot.plot(
        number_of_particles_in_group,
        fraction_of_mass_in_group,
    )
    
    subplot.set_xscale('log')
    subplot.set_yscale('log')
    subplot.set_xlabel('N')
    subplot.set_ylabel('df/d(Log_10 N)')
    
    figure.savefig('x.png')
    pyplot.show()
    
if __name__ in ('__main__', '__plot__'):
    number_of_stars = 10000
    stars = new_cluster(number_of_stars)
    total_mass = stars.mass.sum()
    
    unit_converter = nbody_system.nbody_to_si(total_mass, 1 | units.parsec)
    groups = find_clumps(stars, unit_converter)
    plot_clumps(groups, total_mass)
