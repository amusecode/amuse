"""
Particle colors

This module contains functions to compute colors for stars and SPH particles. These can be used 
for visualisation.
"""

import os.path
import numpy
from amuse.units import units, constants
from amuse.support.exceptions import AmuseException
from amuse.datamodel import ParticlesOverlay


__all__ = ["new_particles_with_color", "new_particles_with_blackbody_color", "mu", "u_from_T", "T_from_u"]


def new_particles_with_color(original_particles, red_function, green_function, blue_function, attributes_names=None):
    """
    Returns new color particles. These are bound to the 'original_particles' in 
    the sense that they share their attributes, but have additional attributes 
    'red', 'green', and 'blue'.
    :argument original_particles: the particles for which the color needs to be computed
    :argument red_function: function that computes red color of a particle
    :argument green_function: function that computes green color of a particle
    :argument blue_function: function that computes blue color of a particle
    """
    original_particles.add_calculated_attribute("red", red_function, attributes_names=attributes_names)
    original_particles.add_calculated_attribute("green", green_function, attributes_names=attributes_names)
    original_particles.add_calculated_attribute("blue", blue_function, attributes_names=attributes_names)
    original_particles.add_vector_attribute("color", ["red", "green", "blue"])
    return original_particles

def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of particles in a gas)
    X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
    x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception("Error in calculating mu: mass fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)

def u_from_T(T, mu=mu(Y=0.25, Z=0.02, x_ion=0.1)):
    """
    Computes internal energy from temperature for a monatomic ideal gas. The default mean
    molecular weight is for solar composition with an ionisation fraction of 0.1.
    """
    return 3.0/2.0 * constants.kB * T / mu

def T_from_u(u, mu=mu(Y=0.25, Z=0.02, x_ion=0.1)):
    """
    Computes temperature from internal energy for a monatomic ideal gas. The default mean
    molecular weight is for solar composition with an ionisation fraction of 0.1.
    """
    return 2.0/3.0 * u * mu / constants.kB


class BlackBodyColorFromTemperature(object):
    def __init__(self):
        self.create_temperature_to_RGB_table()
    
    def input_value_to_temperature_converter(self, temperature):
        return temperature
    
    def red_function(self, input_value):
        temperature = self.input_value_to_temperature_converter(input_value)
        index = numpy.searchsorted(self.temperature.number[1:], temperature.value_in(self.temperature.unit))
        return self.red[index]
    
    def green_function(self, input_value):
        temperature = self.input_value_to_temperature_converter(input_value)
        index = numpy.searchsorted(self.temperature.number[1:], temperature.value_in(self.temperature.unit))
        return self.green[index]
    
    def blue_function(self, input_value):
        temperature = self.input_value_to_temperature_converter(input_value)
        index = numpy.searchsorted(self.temperature.number[1:], temperature.value_in(self.temperature.unit))
        return self.blue[index]
    
    def create_temperature_to_RGB_table(self):
        table_file = os.path.join(os.path.dirname(__file__), 'bbr_color.txt')
        with open(table_file, 'r') as infile:
            temperature, red, green, blue = ([], [], [], [])
            for line in infile.readlines():
                words = line.split()
                if line[0] == "#" or words[2] == "10deg":
                    continue
                temperature.append(float(words[0]))
                red.append(float(words[6]))
                green.append(float(words[7]))
                blue.append(float(words[8]))
        self.red = numpy.array(red)
        self.green = numpy.array(green)
        self.blue = numpy.array(blue)
        self.temperature = temperature | units.K


class BlackBodyColorFromInternalEnergy(BlackBodyColorFromTemperature):
    
    def __init__(self, X=None, Y=0.25, Z=0.02, x_ion=0.1):
        self.create_temperature_to_RGB_table()
        self.mu = mu(X=X, Y=Y, Z=Z, x_ion=x_ion)
    
    def input_value_to_temperature_converter(self, u):
        return T_from_u(u, mu=self.mu)
    

def new_particles_with_blackbody_color(original_particles, **kwargs):
    """
    Returns new color particles. These are bound to the 'original_particles' in 
    the sense that they share their attributes, but have additional attributes 
    'red', 'green', and 'blue'. These colors are based on Mitchell Charity's
    blackbody color datafile (bbr_color.txt, see 
    http://www.vendian.org/mncharity/dir3/blackbody/)
    If the particles have a temperature attribute, the colors are computed from 
    these. Otherwise they will be computed from the gas internal energy, using the 
    T_from_u function, in which case the optional keyword arguments X, Y, 
    Z, and x_ion can be supplied.
    
    :argument original_particles: the particles for which the color needs to be computed
    :argument X: hydrogen abundance for T_from_u converter (default: None, i.e. compute from Y and Z)
    :argument Y: helium abundance for T_from_u converter (default: 0.25)
    :argument Z: metal (everything heavier than helium) abundance for T_from_u converter (default: 0.02)
    :argument x_ion: ionisation fraction for T_from_u converter (default: 0.1)
    """
    if hasattr(original_particles, "temperature"):
        colors = BlackBodyColorFromTemperature(**kwargs)
        attributes_names = ["temperature"]
    elif hasattr(original_particles, "u"):
        colors = BlackBodyColorFromInternalEnergy(**kwargs)
        attributes_names = ["u"]
    else:
        raise AmuseException("The particles need to have 'temperature' or 'u' attributes for deriving black body colors")
    
    return new_particles_with_color(
        original_particles, 
        colors.red_function, colors.green_function, colors.blue_function, 
        attributes_names=attributes_names)

