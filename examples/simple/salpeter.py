import numpy 
from matplotlib import pyplot
from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution
def new_cluster(number_of_stars = 1000):
    masses = new_salpeter_mass_distribution(
        number_of_stars, 
        mass_min = 0.1 | units.MSun,
        mass_max = 125.0 | units.MSun, 
        alpha = -2.35
    )
    nbody_converter = nbody_system.nbody_to_si(masses.sum(), 1 | units.parsec)
    particles = new_plummer_sphere(number_of_stars, nbody_converter)
    particles.mass = masses
    particles.move_to_center()
    return particles

def plot_particles_and_mass_distribution(particles):
    figure = pyplot.figure(figsize= (12,6))
    
    subplot = figure.add_subplot(1, 2, 1)
        
    subplot.scatter(
        particles.x.value_in(units.parsec),
        particles.y.value_in(units.parsec),
        s = particles.mass.value_in(units.MSun),# * len(particles),
        edgecolors = 'red',
        facecolors = 'red'
    )
    
    subplot.set_xlim(-4,4)
    subplot.set_ylim(-4,4)
    subplot.set_xlabel('x (parsec)')
    subplot.set_ylabel('y (parsec)')
    
    subplot = figure.add_subplot(1, 2, 2)
    
    masses = particles.mass.value_in(units.MSun)
    
    bins = 10**numpy.linspace(-1, 2, 100)
    number_of_particles, bin_edges= numpy.histogram(masses, bins = bins)
    
    bin_sizes = bin_edges[1:] - bin_edges[:-1]
    
    y = number_of_particles / bin_sizes
    x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    
    y = y[number_of_particles > 10.0] 
    x = x[number_of_particles > 10.0]
    subplot.scatter(x, y)
    
    c = ((0.1**-1.35) - (125.0**-1.35)) / 1.35
    subplot.plot(x, len(particles)/ c * (x**-2.35))
    
    subplot.set_xscale('log')
    subplot.set_yscale('log')
    
    subplot.set_xlabel(u'M [M\u2299]')
    subplot.set_ylabel('N')
    
    pyplot.show()
    
if __name__ in ('__main__', '__plot__'):
    particles = new_cluster(20000)
    plot_particles_and_mass_distribution(particles)
