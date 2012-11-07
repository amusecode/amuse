"""
Generate a grid of binaries with different, primary mass, mass ratio
and separation and evolve these over time.
"""

from amuse.units import units
from amuse.units import quantities
from amuse import datamodel
from amuse.community.seba.interface import SeBa
from amuse.community.bse.interface import BSE

from matplotlib import pyplot

import numpy

def generate_initial_population_grid(
    min_mass, max_mass, number_of_mass_bins, 
    min_ratio, max_ratio, number_of_ratio_bins,
    min_separation, max_separation, number_of_separation_bins
):
    """
    creates a set of binaries and a set of stars, the grid
    will be divided in equal parts amongst all dimensions (primary mass,
    ratio and separation), the total number of binaries will
    be the product of all bins.
    """
    
    # vector quantity arange takes care of the units, but is otherwhise
    # equal to numpy.arange
    step = (max_mass - min_mass) / number_of_mass_bins
    primary_masses  = quantities.VectorQuantity.arange(min_mass, max_mass+step, step)
    step = (max_ratio - min_ratio) / number_of_ratio_bins
    ratios = numpy.arange(min_ratio, max_ratio + step, step)
    step = (max_separation - min_separation) / number_of_separation_bins
    separations = quantities.VectorQuantity.arange(min_separation, max_separation+step, step)
    
    binaries = datamodel.Particles()
    stars = datamodel.Particles()
    
    for primary_mass in primary_masses:
        for ratio in ratios:
            for separation in separations:
                primary_star = datamodel.Particle()
                primary_star.mass = primary_mass
                
                # we add the particle to the stars set
                # we get a reference to the particle in the stars set
                # back and we want to use that in constructing the
                # binary, so that the binaries all refer to 
                # the stars in the "stars set"
                primary_star = stars.add_particle(primary_star)
                
                secondary_star = datamodel.Particle()
                secondary_star.mass = primary_mass * ratio
                secondary_star = stars.add_particle(secondary_star)
                
                binary = datamodel.Particle()
                binary.eccentricity = 0.0
                binary.semi_major_axis = separation
                binary.child1 = primary_star
                binary.child2 = secondary_star
                
                binaries.add_particle(binary)
                
    return binaries, stars

def evolve_population(binaries, stars, end_time, time_step):
    """ evolve the binary population to end_time and store all
    properties every time_step
    """
    code = SeBa()
    # add the stars first, as the binaries will
    # refer to them
    code.particles.add_particles(stars)
    code.binaries.add_particle(binaries)
    
    time = 0.0 * end_time
    
    code.commit_particles()
    
    channel_from_code_to_model_for_binaries = code.binaries.new_channel_to(binaries)
    channel_from_code_to_model_for_stars = code.particles.new_channel_to(stars)
    
    
    print "start evolving..."
    while time <= end_time:
        time += time_step
        code.evolve_model(time)
        
        print "evolved to time: ", time.as_quantity_in(units.Myr)
        channel_from_code_to_model_for_stars.copy()
        channel_from_code_to_model_for_binaries.copy()
        
        binaries.savepoint(timestamp = time)
        stars.savepoint(timestamp = time)
def make_hr_diagram(stars):
    
    all_luminosities  = [] | units.LSun
    all_temperatures = [] | units.K
    
    for stars_at_time in list(stars.history)[1:]:
        all_luminosities.extend(stars_at_time.luminosity)
        all_temperatures.extend(stars_at_time.temperature)

    pyplot.figure(figsize = (8, 6))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    
    pyplot.scatter(
        all_temperatures.value_in(units.K), 
        all_luminosities.value_in(units.LSun)
    )
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlim(pyplot.xlim()[::-1])
    pyplot.ylim(.01,1.e5)
    pyplot.show()
    

if __name__ == '__main__':
    binaries, stars = generate_initial_population_grid(
        0.1 | units.MSun, 4.0 | units.MSun, 20,
        0.5, 1.0, 2,
        1.0 | units.RSun, 2.0 | units.RSun, 1
    )
    print "generated a population of", len(binaries), "binaries"
    evolve_population(binaries, stars, 1 | units.Gyr, 1 | units.Myr)
    
    make_hr_diagram(stars)
    
