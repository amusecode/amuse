import numpy 
from matplotlib import pyplot
from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
from amuse.community.hermite0.interface import Hermite


def new_cluster(number_of_stars = 1000, radius = None):
    if radius == None:
        radius = (0.5 / number_of_stars) | nbody_system.length
        
    particles = new_plummer_sphere(number_of_stars)
    particles.mass = new_salpeter_mass_distribution_nbody(number_of_stars)
    particles.radius = radius
    particles.move_to_center()
    return particles

def plot_particles_and_highlight_collision(particles, particles1, particles2):
    figure = pyplot.figure()
    subplot = figure.add_subplot(1, 1, 1)
      
    subplot.scatter(
        particles.x.value_in(nbody_system.length),
        particles.y.value_in(nbody_system.length),
        s =  (2 * subplot.bbox.width  * particles.radius.value_in(nbody_system.length)) ** 2 ,
        edgecolors = 'none',
        facecolors = 'red',
        marker = 'o',
    )
    
    subplot.scatter(
        particles1.x.value_in(nbody_system.length),
        particles1.y.value_in(nbody_system.length),
        s = 2 * (2 * subplot.bbox.width  * particles1.radius.value_in(nbody_system.length)) ** 2,
        edgecolors = 'none',
        facecolors = 'black',
        marker = 'o',
    )
    
    subplot.scatter(
        particles2.x.value_in(nbody_system.length),
        particles2.y.value_in(nbody_system.length),
        s = 2 * (2 * subplot.bbox.width  *  particles2.radius.value_in(nbody_system.length)) ** 2,
        edgecolors = 'none',
        facecolors = 'black',
        marker = 'o',
    )
    
    subplot.set_xlim(-1,1)
    subplot.set_ylim(-1,1)
    subplot.set_xlabel('x (nbody length)')
    subplot.set_ylabel('y (nbody length)')
    
    pyplot.show()
    
if __name__ in ('__main__', '__plot__'):
    numpy.random.seed(1212)
    
    particles = new_cluster(128)
    code = Hermite()
    code.particles.add_particles(particles)
    stopping_condition = code.stopping_conditions.collision_detection
    stopping_condition.enable()
    
    code.evolve_model(4 | nbody_system.time)
    if not stopping_condition.is_set():
        raise Exception("No stopping collision detected in the given timeframe.")
        
    plot_particles_and_highlight_collision(particles, stopping_condition.particles(0), stopping_condition.particles(1) )

