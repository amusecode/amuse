import math
import numpy

from amuse.datamodel import rotation
from amuse.datamodel import Particles
from amuse.units import nbody_system

def new_binary_from_orbital_elements(
        mass1,
        mass2,
        semi_major_axis, 
        eccentricity = 0,
        true_anomaly = 0, 
        inclination = 0,
        longitude_of_the_ascending_node = 0,
        argument_of_periapsis = 0,
        G = nbody_system.G
    ):
    inclination = numpy.radians(inclination)
    argument_of_periapsis = numpy.radians(argument_of_periapsis)
    longitude_of_the_ascending_node = numpy.radians(longitude_of_the_ascending_node)
    true_anomaly = numpy.radians(true_anomaly)

    q   = mass2 / mass1
    mu  = G * (mass2 + mass1)
    
    semilatus_rectum  = semi_major_axis*(1-eccentricity**2)
    radius            = semilatus_rectum/(1+eccentricity*numpy.cos(true_anomaly))
    velocity          = numpy.sqrt(mu/semilatus_rectum)*eccentricity*numpy.sin(true_anomaly)
    
    specific_relative_angular_momentum = numpy.sqrt(mu*semilatus_rectum)
    
    r = radius.new_zeros_array(3)
    v = velocity.new_zeros_array(3)
    
    r[0] = radius * numpy.cos(true_anomaly)
    r[1] = radius * numpy.sin(true_anomaly) 
    
    v[0] = velocity * numpy.cos(true_anomaly) - specific_relative_angular_momentum / radius * numpy.sin(true_anomaly)
    v[1] = velocity * numpy.sin(true_anomaly) + specific_relative_angular_momentum / radius * numpy.cos(true_anomaly)

    result = Particles(2)
    result.position = radius.new_zeros_array(3)
    result.velocity = velocity.new_zeros_array(3)
    result[0].mass = mass1
    result[1].mass = mass2

    
    result[1].position = r
    result[1].velocity = v
    rotation.rotate(result, argument_of_periapsis, inclination, longitude_of_the_ascending_node)
    return result
