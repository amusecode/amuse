"""
Generates a grid of binaries with different, primary mass, mass ratio
and separation and evolves these over time.
"""

from amuse.units import units
from amuse.units import quantities
from amuse import datamodel
from amuse.community.seba.interface import SeBa
from amuse.community.bse.interface import BSE

from matplotlib import pyplot

import numpy

def create_binary(stars, binaries, primary_mass, mass_ratio, separation, eccentricity):
    primary_star = datamodel.Particle()
    primary_star.mass = primary_mass
    
    # we add the particle to the stars set
    # and we get a reference to the particle in the stars set
    # back
    # we want to use that star in constructing the
    # binary, so that the binaries all refer to 
    # the stars in the "stars set"
    primary_star = stars.add_particle(primary_star)
    
    secondary_star = datamodel.Particle()
    secondary_star.mass = primary_mass * mass_ratio
    secondary_star = stars.add_particle(secondary_star)
    
    binary = datamodel.Particle()
    binary.eccentricity = eccentricity
    binary.semi_major_axis = separation
    binary.child1 = primary_star
    binary.child2 = secondary_star
    
    binaries.add_particle(binary)
def generate_initial_population_grid(
    min_mass, max_mass, number_of_mass_bins, 
    min_ratio, max_ratio, number_of_ratio_bins,
    min_separation, max_separation, number_of_separation_bins,
    min_eccentricity, max_eccentricity, number_of_eccentricity_bins
):
    """
    creates a set of binaries and a set of stars, the grid
    will be divided in equal parts amongst all dimensions (primary mass,
    ratio and separation), the total number of binaries will
    be the product of all bins.
    """
    
    # vector quantity arange takes care of the units, but is otherwhise
    # equal to numpy.arange
    primary_masses  = quantities.linspace(min_mass, max_mass, number_of_mass_bins)
    mass_ratios = quantities.linspace(min_ratio, max_ratio, number_of_ratio_bins)
    separations = quantities.linspace(min_separation, max_separation, number_of_separation_bins)
    eccentricities = quantities.linspace(min_eccentricity, max_eccentricity, number_of_eccentricity_bins)
    
    binaries = datamodel.Particles()
    stars = datamodel.Particles()
    
    for primary_mass in primary_masses:
        for mass_ratio in mass_ratios:
            for separation in separations:
                for eccentricity in eccentricities :
                    create_binary(
                        stars, 
                        binaries, 
                        primary_mass, 
                        mass_ratio, 
                        separation, 
                        eccentricity
                    )
                
    return binaries, stars

def evolve_population(binaries, stars, end_time, time_step):
    code = SeBa()
    
    # add the stars first, as the binaries will
    # refer to them
    code.particles.add_particles(stars)
    code.binaries.add_particle(binaries)
    
    
    channel_from_code_to_model_for_binaries = code.binaries.new_channel_to(binaries)
    channel_from_code_to_model_for_stars = code.particles.new_channel_to(stars)
    
    #we evolve in steps of timestep, just to get some feedback
    print "start evolving..."
    time = 0.0 * end_time
    while time < end_time:
        time += time_step
        code.evolve_model(time)
        print "evolved to time: ", time.as_quantity_in(units.Myr)
        
    channel_from_code_to_model_for_stars.copy()
    channel_from_code_to_model_for_binaries.copy()
    
def make_hr_diagram(binaries):
    pyplot.figure(figsize = (8, 6))
    pyplot.title('Binary population', fontsize=12)
    separation = binaries.semi_major_axis.value_in(units.RSun)
    eccentricity = binaries.eccentricity
    pyplot.hexbin(
        separation,
        eccentricity,
        gridsize = 25,
        bins = 'log',
        extent = (0.0, 2.0, 0.0, 0.9)
    )
    pyplot.show()
    

if __name__ == '__main__':
    print "generating a binary population..."
    
    binaries, stars = generate_initial_population_grid(
        5 | units.MSun, 10.0 | units.MSun, 1, #mass range
        0.5, 0.9, 1,  #mass ratios range
        0.2 | units.RSun, 2.0 | units.RSun, 60, #semi major axis range
        0.0 , 1.0, 60 #eccentricity range
    )
    
    print "generated a population of", len(binaries), "binaries"
    
    evolve_population(binaries, stars,  1 | units.Gyr, 500 | units.Myr)
    
    make_hr_diagram(binaries)
    
