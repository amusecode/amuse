"""
Salpeter initial mass function (IMF)

This module contains functions used to generate realisations of the 
Salpeter IMF.
"""

import numpy
import os
import sys

from amuse.units import nbody_system
from amuse.units import units

__all__ = ["new_salpeter_mass_distribution", "new_salpeter_mass_distribution_nbody"]

class SalpeterIMF(object):
    def __init__(self, mass_min = 0.1 | units.MSun, mass_max = 125 | units.MSun, alpha = -2.35):
        self.mass_min = mass_min
        self.mass_max = mass_max
        self.alpha = alpha
        self.random = numpy.random
    
    def mass_mean(self):
        alpha1 = self.alpha + 1
        alpha2 = self.alpha + 2
        l1 = pow(self.mass_min, alpha1)
        l2 = pow(self.mass_min, alpha2)
        u1 = pow(self.mass_max, alpha1)
        u2 = pow(self.mass_max, alpha2)
        return ((u2 - l2) * alpha1) / ((u1 - l1) * alpha2)
        
    def mass(self, random_number):
        alpha1 = self.alpha + 1
        factor = (pow(self.mass_max / self.mass_min, alpha1) - 1.0 )
        return self.mass_min * (pow(1.0 + (factor * random_number), 1.0 / alpha1))
        
    def next_mass(self,N=1):
        return self.mass(self.random.random(N))
    
    def next_set(self, number_of_stars):
        set_of_masses=self.next_mass(number_of_stars).in_(self.mass_min.unit)
        total_mass=set_of_masses.sum()
        return (total_mass, set_of_masses)
    

def new_salpeter_mass_distribution(number_of_particles, *list_arguments, **keyword_arguments):
    """Returns a salpeter mass distribution in SI units.
    
    :argument alpha: the dimensionless exponent of the Salpeter function (defaults to -2.35)
    """
    uc = SalpeterIMF(*list_arguments, **keyword_arguments)
    return uc.next_mass(number_of_particles).as_quantity_in(uc.mass_min.unit)
    
def new_salpeter_mass_distribution_nbody(number_of_particles, **keyword_arguments):
    """Returns a salpeter mass distribution in nbody masses.
    All masses will be scaled so that the total mass is always 1.0.
    
    :argument alpha: the dimensionless exponent of the Salpeter function (defaults to -2.35)
    """
    if not 'mass_min' in keyword_arguments:
        keyword_arguments['mass_min'] = 0.1 | nbody_system.mass
        
    if not 'mass_max' in keyword_arguments:
        keyword_arguments['mass_max'] = 125 | nbody_system.mass
        
    uc = SalpeterIMF(**keyword_arguments)
    total_mass, result = uc.next_set(number_of_particles)
    result *=  (1.0 | total_mass.unit) / total_mass
    return result
    
