import sys
import numpy
import os

from amuse.support.units import units

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
        one = 1.0 | units.none
        alpha1 = self.alpha + 1
        factor = (pow(self.mass_max / self.mass_min , alpha1) - one )
        return self.mass_min * (pow(one + (factor * random_number), 1.0 / alpha1))
        
    def next_mass(self):
        return self.mass(self.random.random())
        
    def next_set(self, number_of_stars):
        set_of_masses = self.mass_min.unit.new_quantity(numpy.zeros(number_of_stars))
        total_mass = self.mass_min.unit.new_quantity(0.0)
        for i in range(number_of_stars):
           mass = self.next_mass()
           set_of_masses[i] = mass
           total_mass += mass
        return (total_mass, set_of_masses)
        
