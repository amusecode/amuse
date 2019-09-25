import os
import os.path
import shutil
import numpy

from amuse.lab import *
from amuse.community.mesa.interface import MESA as stellar_evolution_code
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model

from matplotlib import pyplot

def plot_clumps(groups):
    number_of_particles_in_group = []
    fraction_of_mass_in_group =  []

#        number_of_particles_in_group.append(len(group))
#        fraction = (group.mass.sum()/total_mass)
#        fraction_of_mass_in_group.append(fraction)
    
    print "N=", len(groups)
    ci = ['r', 'b', 'g', 'k']
    figure = pyplot.figure(figsize= (12,6))
    i = 0
    alpha = 1
    sizes = 50
    for group in groups:
        pyplot.scatter(group.x.value_in(units.RSun), group.y.value_in(units.RSun), sizes, ci[i], edgecolors = "none", alpha = alpha)


#        pyplot.scatter(
#            group.x.value_in(units.RSun),
#            group.y.value_in(units.RSun),
#            s = 1,#group.mass.value_in(units.MSun),
#            c = ci[i]
#            )
        i+=1
    pyplot.xlabel('x (AU)')
    pyplot.ylabel('y (A*)')
#    pyplot.xlim(-30, 30)
#   pyplot.ylim(-30, 30)
    pyplot.show()

def find_clumps(particles, unit_converter):

    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()

    mean_densty = hop.particles.density.mean() 
    hop.parameters.peak_density_threshold = mean_densty
    hop.parameters.saddle_density_threshold = 0.99*mean_densty
    hop.parameters.outer_density_threshold = 0.01*mean_densty
#    print "Peak density treshold:", 

    hop.do_hop()
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    hop.stop()
    return result


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", 
                      default = "hydro_triple_gas.hdf5",
                      help="input filename [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="time", type="float", 
                      default = 1|units.Myr,
                      help="evolution time scale [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float", 
                      default = 3|units.MSun,
                      help="stellar mass [%default]")
    return result

if __name__ == "__main__":
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments = new_option_parser().parse_args()

    print "initialize star"
    stellar_evolution = EVtwin()
    stellar_evolution.particles.add_particle(Particle(mass=o.Mprim))
    stellar_evolution.evolve_model(o.time)
    particles = convert_stellar_model_to_SPH(
        stellar_evolution.particles[0], 
        500, 
        seed=12345
    ).gas_particles
    stellar_evolution.stop()
    from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits

    print "convert star to SPH"
    converter = nbody_system.nbody_to_si(1|units.hour, 1|units.RSun)
    hydrodynamics = Gadget2(converter)
    hydrodynamics.gas_particles.add_particles(particles)
    hydrodynamics.evolve_model(1.0|units.s)
    hydrodynamics.gas_particles.copy_values_of_attributes_to(["density", "u", "pressure"], particles)
    hydrodynamics.stop()

    print "convert SPH to stellar model"
    model = convert_SPH_to_stellar_model(particles)

    #stellar_evolution = MESA()
    print "initiate stellar evolution model"
    stellar_evolution = EVtwin(redirect="none")    
    stellar_evolution.new_particle_from_model(model, 0.0|units.Myr)
    print "star:", stellar_evolution.particles
    print "evolve star"
    stellar_evolution.evolve_model(2*o.time)
    print stellar_evolution.particles

    



