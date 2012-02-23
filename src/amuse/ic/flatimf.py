"""
Flat initial mass function (IMF)

This module contains functions used to generate realisations of a 
logarithmically flat IMF.
"""

import numpy
from amuse.units import units, nbody_system

__all__ = ["new_flat_mass_distribution", "new_flat_mass_distribution_nbody"]

class FlatIMF(object):
    def __init__(self, mass_min = 0.1 | units.MSun, mass_max = 125.0 | units.MSun):
        self.mass_min = mass_min
        self.mass_max = mass_max
        self.random = numpy.random
    
    def mass_mean(self):
        return ((self.mass_max - self.mass_min) / 
            numpy.log(self.mass_max / self.mass_min)).as_quantity_in(self.mass_min.unit)
        
    def mass(self, random_number):
        return self.mass_min * (self.mass_max / self.mass_min)**random_number
        
    def next_mass(self,N=1):
        return self.mass(self.random.random(N))
    
    def next_set(self, number_of_stars):
        set_of_masses=self.next_mass(number_of_stars).in_(self.mass_min.unit)
        total_mass=set_of_masses.sum()
        return (total_mass, set_of_masses)
    

def new_flat_mass_distribution(number_of_particles, *list_arguments, **keyword_arguments):
    """Returns a logarithmically flat mass distribution in SI units.
    
    :argument mass_min: the minimum mass of the distribution (defaults to 0.1 MSun)
    :argument mass_max: the maximum mass of the distribution (defaults to 125.0 MSun)
    """
    uc = FlatIMF(*list_arguments, **keyword_arguments)
    return uc.next_mass(number_of_particles).as_quantity_in(uc.mass_min.unit)
    
def new_flat_mass_distribution_nbody(number_of_particles, **keyword_arguments):
    """Returns a logarithmically flat mass distribution in nbody masses.
    All masses will be scaled so that the total mass is always 1.0.
    
    :argument mass_min: the minimum mass of the distribution (defaults to 0.1)
    :argument mass_max: the maximum mass of the distribution (defaults to 125.0)
    """
    if not 'mass_min' in keyword_arguments:
        keyword_arguments['mass_min'] = 0.1 | nbody_system.mass
        
    if not 'mass_max' in keyword_arguments:
        keyword_arguments['mass_max'] = 125.0 | nbody_system.mass
        
    uc = FlatIMF(**keyword_arguments)
    total_mass, result = uc.next_set(number_of_particles)
    result *=  (1.0 | total_mass.unit) / total_mass
    return result
    
