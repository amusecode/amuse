"""
Plummer model generator

This module contains a function used to create Plummer (1911) models, which 
follow a spherically symmetric density profile of the form:
rho = c * (1 + r**2)**(-5/2)
"""

import numpy
import numpy.random
import random

from math import pi
from math import sqrt

from amuse.units import nbody_system
from amuse import datamodel

__all__ = ["new_plummer_sphere"]

class MakePlummerModel(object):
    def __init__(self, number_of_particles, convert_nbody = None, radius_cutoff = 22.8042468, mass_cutoff = 0.999,
            do_scale = False, random_state = None):
        self.number_of_particles = number_of_particles
        self.convert_nbody = convert_nbody
        self.mass_cutoff = min(mass_cutoff, self.calculate_mass_cuttof_from_radius_cutoff(radius_cutoff))
        self.do_scale = do_scale
        if not random_state == None:
            print "DO NOT USE RANDOM STATE"
        self.random_state = None

    def calculate_mass_cuttof_from_radius_cutoff(self, radius_cutoff):
        if radius_cutoff > 99999:
            return 0.999
        scale_factor = 16.0 / (3.0 * pi)
        rfrac = radius_cutoff * scale_factor
        denominator = pow(1.0 + rfrac ** 2, 1.5)
        numerator = rfrac ** 3
        return numerator/denominator

    def calculate_radius(self, index):
        mass_min = (index * self.mass_cutoff) / self.number_of_particles
        mass_max = ((index+1) * self.mass_cutoff) / self.number_of_particles
        random_mass_fraction = numpy.random.uniform(mass_min, mass_max)
        radius = 1.0 / sqrt( pow (random_mass_fraction, -2.0/3.0) - 1.0)
        return radius

    def calculate_radius_uniform_distribution(self):
        return 1.0 /  numpy.sqrt( numpy.power(numpy.random.uniform(0,self.mass_cutoff,(self.number_of_particles,1)), -2.0/3.0) - 1.0)

    def new_positions_spherical_coordinates(self):
        pi2 = pi * 2
        result = numpy.zeros((self.number_of_particles, 3))
        radius = self.calculate_radius_uniform_distribution()
        theta = numpy.arccos(numpy.random.uniform(-1.0,1.0, (self.number_of_particles,1)))
        phi = numpy.random.uniform(0.0,pi2, (self.number_of_particles,1))
        return (radius,theta,phi)

    def new_velocities_spherical_coordinates(self, radius):
        pi2 = pi * 2
        x,y = self.new_xy_for_velocity()
        velocity = x * sqrt(2.0) * numpy.power( 1.0 + radius*radius, -0.25)
        theta = numpy.arccos(numpy.random.uniform(-1.0,1.0, (self.number_of_particles,1)))
        phi = numpy.random.uniform(0.0,pi2, (self.number_of_particles,1))
        return (velocity,theta,phi)

    def coordinates_from_spherical(self, radius, theta, phi):
        x = radius * numpy.sin( theta ) * numpy.cos( phi )
        y = radius * numpy.sin( theta ) * numpy.sin( phi )
        z = radius * numpy.cos( theta )
        return (x,y,z)

    def new_xy_for_velocity(self):
        number_of_selected_items = 0
        selected_values_for_x = numpy.zeros(0)
        selected_values_for_y = numpy.zeros(0)
        while (number_of_selected_items < self.number_of_particles):
            x = numpy.random.uniform(0,1.0, (self.number_of_particles-number_of_selected_items))
            y = numpy.random.uniform(0,0.1, (self.number_of_particles-number_of_selected_items))
            g = (x**2) * numpy.power(1.0 - x**2, 3.5)
            compare = y <= g
            selected_values_for_x = numpy.concatenate((selected_values_for_x, x.compress(compare)))
            selected_values_for_y= numpy.concatenate((selected_values_for_x, y.compress(compare)))
            number_of_selected_items = len(selected_values_for_x)
        return numpy.atleast_2d(selected_values_for_x).transpose(), numpy.atleast_2d(selected_values_for_y).transpose()

    def new_model(self):
        m = numpy.zeros((self.number_of_particles,1)) + (1.0 / self.number_of_particles)
        radius, theta, phi = self.new_positions_spherical_coordinates()
        position =  numpy.hstack(self.coordinates_from_spherical(radius, theta, phi))
        radius, theta, phi = self.new_velocities_spherical_coordinates(radius)
        velocity = numpy.hstack(self.coordinates_from_spherical(radius, theta, phi))
        position = position / 1.695
        velocity = velocity / sqrt(1 / 1.695)
        return (m, position, velocity)

    @property
    def result(self):
        masses, positions, velocities = self.new_model()
        result = datamodel.Particles(self.number_of_particles)
        result.mass = nbody_system.mass.new_quantity(numpy.hstack(masses))
        result.position = nbody_system.length.new_quantity(positions)
        result.velocity = nbody_system.speed.new_quantity( velocities)
        result.radius = 0 | nbody_system.length

        result.move_to_center()
        if self.do_scale:
            result.scale_to_standard()

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()

        return result

def new_plummer_sphere(number_of_particles, *list_arguments, **keyword_arguments):
    """
    Create a plummer sphere with the given number of particles. Returns
    a set of stars with equal mass and positions and velocities distributed
    to fit a plummer star distribution model. The model is centered around the
    origin. Positions and velocities are optionally scaled such that the kinetic and
    potential energies are 0.25 and -0.5 in nbody-units, respectively.

    :argument number_of_particles: Number of particles to include in the plummer sphere
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument radius_cutoff: Cutoff value for the radius (defaults to 22.8042468)
    :argument mass_cutoff: Mass percentage inside radius of 1
    :argument do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)
    """
    uc = MakePlummerModel(number_of_particles, *list_arguments, **keyword_arguments)
    return uc.result
