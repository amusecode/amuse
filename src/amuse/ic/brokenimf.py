"""
Broken power-law initial mass function (IMF)

This module contains functions used to generate realisations of multiple-part 
power-law IMFs such as Kroupa (2001), Scalo (1986), Miller & Scalo (1979).
"""

import numpy
import numpy.random

from amuse.units import units
from amuse.units.quantities import AdaptingVectorQuantity

__all__ = ["MultiplePartIMF", "new_broken_power_law_mass_distribution", "new_scalo_mass_distribution", 
    "new_miller_scalo_mass_distribution", "new_kroupa_mass_distribution"]

class MultiplePartIMF(object):
    def __init__(self, mass_boundaries = [0.1, 125.0] | units.MSun, mass_max = None, alphas = [-2.35], random=True):
        self.mass_boundaries = mass_boundaries
        if not mass_max is None:
            self.mass_boundaries[-1] = mass_max
        self.alphas = numpy.array(alphas)
        
        if random:
            self.random = numpy.random.random
        else:
            self.random = self.evenly_distributed
        
        self.number_of_bins = len(self.alphas)
        self.fraction_per_bin = self.calculate_fraction_per_bin()
        self.cumulative_fractions = numpy.array([sum(self.fraction_per_bin[:i]) for i in range(self.number_of_bins+1)])
        self.cumulative_fractions[-1] = 1.0 # In case of round-off errors
        self.factors = pow(self.mass_boundaries[1:] / self.mass_boundaries[:-1], self.alphas + 1.0) - 1.0
        self.inv_alpha1s = numpy.array([numpy.inf if alpha==-1 else (1.0 / (alpha + 1.0)) for alpha in self.alphas])
    
    def evenly_distributed(self, size=1):
        return numpy.linspace(0.0, 1.0, size)
    
    def calculate_fraction_per_bin(self):
        result = []
        mass_boundaries = self.mass_boundaries.value_in(units.MSun)
        for i, alpha in enumerate(self.alphas):
            if alpha == -1:
                factor = numpy.log(mass_boundaries[i+1] / mass_boundaries[i])
            else:
                factor = (mass_boundaries[i+1]**(alpha+1) - mass_boundaries[i]**(alpha+1)) / (alpha+1)
            j=0
            for j in range(self.number_of_bins - i - 1):
                factor *= mass_boundaries[-j-2]**(self.alphas[-j-1]-self.alphas[-j-2])
            
            result.append(factor)
        total = sum(result, 0.0)
        return numpy.array(result)/total
    
    def mass_mean(self):
        result = 0 | units.MSun
        for i in range(self.number_of_bins):
            a1 = self.alphas[i] + 1
            a2 = self.alphas[i] + 2
            m_low = self.mass_boundaries[i]
            m_high = self.mass_boundaries[i+1]
            
            if abs(a1) < 1.0e-10:
                mean_of_bin = (m_high**a2 - m_low**a2) / (a2 * numpy.log(m_high/m_low))
            elif abs(a2) < 1.0e-10:
                mean_of_bin = (a1 / (m_high**a1 - m_low**a1)) * numpy.log(m_high/m_low)
            else:
                mean_of_bin = (a1 / (m_high**a1 - m_low**a1)) * ((m_high**a2 - m_low**a2) / a2)
            result += self.fraction_per_bin[i] * mean_of_bin
        return result
        
    def mass(self, random_numbers):
        indices = numpy.searchsorted(self.cumulative_fractions[:-1], random_numbers, 'right') - 1
        scaled = ((random_numbers - self.cumulative_fractions[indices]) / 
            (self.cumulative_fractions[indices+1] - self.cumulative_fractions[indices]))
        
        result = numpy.empty_like(random_numbers)
        zerodiv = self.alphas[indices]==-1
        normal = numpy.logical_not(zerodiv)
        result[zerodiv] = pow(self.mass_boundaries[1:][indices[zerodiv]] / self.mass_boundaries[:-1][indices[zerodiv]], scaled[zerodiv])
        result[normal] = pow(1.0 + self.factors[indices[normal]] * scaled[normal], self.inv_alpha1s[indices[normal]])
        return self.mass_boundaries[:-1][indices] * result
        
    def next_mass(self, number_of_stars=1):
        return self.mass(self.random(number_of_stars)).as_quantity_in(self.mass_boundaries.unit)
    
    def next_set(self, number_of_stars):
        set_of_masses = self.next_mass(number_of_stars)
        total_mass = set_of_masses.sum()
        return (total_mass, set_of_masses)
    

def new_broken_power_law_mass_distribution(number_of_particles, *list_arguments, **keyword_arguments):
    """Returns a broken (multiple-part) power-law mass distribution in SI units.
    
    :argument mass_boundaries: the boundaries of the mass ranges (default: [0.1, 125.0] MSun)
    :argument alphas: the exponents of each mass range (default: [-2.35])
    """
    uc = MultiplePartIMF(*list_arguments, **keyword_arguments)
    return uc.next_mass(number_of_particles)

def new_scalo_mass_distribution(number_of_particles, mass_max=None, random=True):
    """Returns a Scalo (1986) mass distribution in SI units, with mass ranges:
        [0.10, 0.18, 0.42, 0.62, 1.18, 3.5, 125.0] MSun,
    and power-law exponents of each mass range:
        [1.6, -1.01, -2.75, -2.08, -3.5, -2.63]
    
    :argument mass_max: the cut-off mass (defaults to 125.0 MSun)
    """
    uc = MultiplePartIMF(
        mass_boundaries = [0.10, 0.18, 0.42, 0.62, 1.18, 3.5, 125.0] | units.MSun, 
        mass_max = mass_max, alphas = [1.6, -1.01, -2.75, -2.08, -3.5, -2.63], random=random)
    return uc.next_mass(number_of_particles)

def new_miller_scalo_mass_distribution(number_of_particles, mass_max=None, random=True):
    """Returns a Miller & Scalo (1979) mass distribution in SI units, with mass ranges:
        [0.1, 1.0, 2.0, 10.0, 125.0] MSun,
    and power-law exponents of each mass range:
        [-1.25, -2.0, -2.3, -3.3]
    
    :argument mass_max: the cut-off mass (defaults to 125.0 MSun)
    """
    uc = MultiplePartIMF(
        mass_boundaries = [0.1, 1.0, 2.0, 10.0, 125.0] | units.MSun, 
        mass_max = mass_max, alphas = [-1.25, -2.0, -2.3, -3.3], random=random)
    return uc.next_mass(number_of_particles)

def new_kroupa_mass_distribution(number_of_particles, mass_max=None, random=True):
    """Returns a Kroupa (2001) mass distribution in SI units, with mass ranges:
        [0.01, 0.08, 0.5, 100.0] MSun,
    and power-law exponents of each mass range:
        [-0.3, -1.3, -2.3]
    
    :argument mass_max: the cut-off mass (defaults to 100.0 MSun)
    """
    uc = MultiplePartIMF(
        mass_boundaries = [0.01, 0.08, 0.5, 100.0] | units.MSun, 
        mass_max = mass_max, alphas = [-0.3, -1.3, -2.3], random=random)
    return uc.next_mass(number_of_particles)


