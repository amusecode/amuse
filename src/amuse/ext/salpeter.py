import sys
import numpy
import random
import os

from amuse.support.units import units

class SalpeterIMF(object):
    def __init__(self, mass_min = 0.1 | units.MSun, mass_max = 125 | units.MSun, alpha = -2.35):
        self.mass_min = mass_min.as_quantity_in(units.MSun)
        self.mass_max = mass_max.as_quantity_in(units.MSun)
        self.alpha = alpha
        self.random = random.Random()
        self.random.seed()
    
    def mass_mean(self):
        alpha1 = self.alpha + 1
        alpha2 = self.alpha + 2
        l1 = pow(self.mass_min.value_in(units.MSun), alpha1)
        l2 = pow(self.mass_min.value_in(units.MSun), alpha2)
        u1 = pow(self.mass_max.value_in(units.MSun), alpha1)
        u2 = pow(self.mass_max.value_in(units.MSun), alpha2)
        return ((u2 - l2) * alpha1) / ((u1 - l1) * alpha2) | units.MSun
        
    def mass(self, random_number):
        alpha1 = self.alpha + 1
        factor = (pow(self.mass_max.value_in(units.MSun) / self.mass_min.value_in(units.MSun) , alpha1) - 1.0)
        return self.mass_min.value_in(units.MSun) * (pow(1 + (factor * random_number), 1.0 / alpha1)) | units.MSun
        
    def next_mass(self):
        return self.mass(self.random.random())
        
    def next_set(self, number_of_stars):
        set_of_masses = numpy.zeros(number_of_stars)
        total_mass = 0.0 | units.MSun
        for i in range(number_of_stars):
           mass = self.next_mass()
           set_of_masses[i] = mass.value_in(units.MSun)
           total_mass += mass
        return (total_mass, units.MSun.new_quantity(set_of_masses))
        
