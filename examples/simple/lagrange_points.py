"""
Make a contour plot of the effective potential of the Sun-Earth system (left) 
and a system with a 10000 times more massive Earth with the Lagrangian points visible.
"""

import numpy
from matplotlib import pyplot
from amuse.plot import xlabel, ylabel, effective_iso_potential_plot

from amuse.units import units, constants, nbody_system
from amuse.community.hermite0.interface import Hermite
from amuse.datamodel import Particles

def new_sun_earth_system():
    particles = Particles(2)
    particles.mass = [1, 3.0024584e-6] | units.MSun
    particles.position = [[0, 0, 0], [1.0, 0, 0]] | units.AU
    particles.velocity = [0, 0, 0] | units.km / units.s
    particles[1].vy = (constants.G * particles.total_mass() / particles[1].x).sqrt()
    return particles

def setup_gravity_code():
    converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.AU)
    gravity = Hermite(converter)
    gravity.parameters.epsilon_squared = 1e-6 | units.RSun**2
    gravity.particles.add_particles(new_sun_earth_system())
    return gravity

def make_effective_iso_potential_plot(gravity_code):
    omega = (constants.G * gravity_code.particles.total_mass() / (1.0|units.AU**3)).sqrt()
    center_of_mass = gravity_code.particles.center_of_mass()[:2]
    figure = pyplot.figure(figsize = (9, 8))
    current_axes = pyplot.subplot(2, 2, 1)
    current_axes.set_aspect("equal", adjustable = "box")
    effective_iso_potential_plot(gravity_code, 
        omega, 
        center_of_rotation = center_of_mass,
        fraction_screen_filled=0.7)
    xlabel('x')
    ylabel('y')
    
    current_axes = pyplot.subplot(2, 2, 3)
    current_axes.set_aspect("equal", adjustable = "box")
    effective_iso_potential_plot(gravity_code, 
        omega, 
        center_of_rotation = center_of_mass,
        xlim=[0.9, 1.1]|units.AU, 
        ylim=[-0.1, 0.1]|units.AU,
        number_of_contours=20)
    xlabel('x')
    ylabel('y')
    
    gravity_code.particles[1].mass *= 10000
    omega = (constants.G * gravity_code.particles.total_mass() / (1.0|units.AU**3)).sqrt()
    center_of_mass = gravity_code.particles.center_of_mass()[:2]
    current_axes = pyplot.subplot(2, 2, 2)
    current_axes.set_aspect("equal", adjustable = "box")
    effective_iso_potential_plot(gravity_code, 
        omega, 
        center_of_rotation = center_of_mass,
        number_of_contours=20,
        fraction_screen_filled=0.7)
    xlabel('x')
    ylabel('y')
    
    current_axes = pyplot.subplot(2, 2, 4)
    current_axes.set_aspect("equal", adjustable = "box")
    effective_iso_potential_plot(gravity_code, 
        omega, 
        center_of_rotation = center_of_mass,
        xlim=[0.6, 1.4]|units.AU, 
        ylim=[-0.4, 0.4]|units.AU,
        number_of_contours=20,
        fraction_screen_filled=0.9)
    xlabel('x')
    ylabel('y')
    pyplot.show()   
    
if __name__ in ('__main__', '__plot__'):
    gravity = setup_gravity_code()
    make_effective_iso_potential_plot(gravity)
    gravity.stop()
