import os
import os.path
import shutil
import numpy

from amuse.lab import *
from amuse.community.mesa.interface import MESA as stellar_evolution_code
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model

from prepare_figure import single_frame
from distinct_colours import get_distinct

from matplotlib import pyplot

def plot_clumps(groups):
    number_of_particles_in_group = []
    fraction_of_mass_in_group =  []

#        number_of_particles_in_group.append(len(group))
#        fraction = (group.mass.sum()/total_mass)
#        fraction_of_mass_in_group.append(fraction)
    
    print "N=", len(groups)
    ci = ['r', 'b', 'g', 'k']
    figure = pyplot.figure(figsize=(12,6))
    i = 0
    alpha = 1
    sizes = 50
    for group in groups:
        pyplot.scatter(group.x.value_in(units.RSun),
                       group.y.value_in(units.RSun),
                       sizes, ci[i], edgecolors = "none", alpha = alpha)
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
    # print "Peak density threshold:", 

    hop.do_hop()
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    hop.stop()
    return result

def evolve_single_star(mass, tend):
        star = Particles(1)
        star.mass = mass
        stellar_evolution = MESA()
        stellar_evolution.particles.add_particles(star)
        time = [] | units.Myr
        mass = [] | units.MSun
        radius = [] | units.RSun
        temperature = [] | units.K
        luminosity = [] | units.LSun
        stellar_type = [] 
        while stellar_evolution.model_time<tend:
            stellar_evolution.evolve_model()
            time.append(stellar_evolution.model_time)
            mass.append(stellar_evolution.particles[0].mass)
            radius.append(stellar_evolution.particles[0].radius)
            temperature.append(stellar_evolution.particles[0].temperature)
            luminosity.append(stellar_evolution.particles[0].luminosity)
            stellar_type.append(stellar_evolution.particles[0].stellar_type)
            print "Time=", time[-1], stellar_type[-1], mass[-1], radius[-1], \
                  temperature[-1].in_(units.K), luminosity[-1].in_(units.LSun)
            if stellar_type[-1] >= 4 | units.stellar_type:
                break

        stellar_evolution.stop()
        return time, stellar_type, mass, radius, temperature, luminosity

def merge_two_stars_sph_and_evolve(Mprim, Msec, tcoll, tend):
    stars = Particles(2)
    stars[0].mass = Mprim
    stars[1].mass = Msec
    stellar = EVtwin()
    stellar.particles.add_particle(stars[0])
    stellar.particles.add_particle(stars[1])
    
    time = [] | units.Myr
    mass = [] | units.MSun
    radius = [] | units.RSun
    temperature = [] | units.K
    luminosity = [] | units.LSun
    while stellar.model_time < tcoll:
        stellar.evolve_model()
        time.append(stellar.model_time)
        mass.append(stellar.particles[0].mass)
        radius.append(stellar.particles[0].radius)
        temperature.append(stellar.particles[0].temperature)
        luminosity.append(stellar.particles[0].luminosity)
        print "Time=", time[-1], mass[-1], radius[-1], \
              temperature[-1].in_(units.K), luminosity[-1].in_(units.LSun)
    n_normal = len(time)

    print stars
    Nprim = int(100*stellar.particles[0].mass.value_in(units.MSun))
    mgas = stellar.particles[0].mass/Nprim
    Nsec = int(stellar.particles[1].mass/mgas)
    print "N gas=", Nprim, Nsec
    sph_primary = convert_stellar_model_to_SPH(
        stellar.particles[0],
        Nprim, 
        seed=12345
    ).gas_particles
    sph_secondary = convert_stellar_model_to_SPH(
        stellar.particles[0], 
        Nsec, 
        seed=12345
    ).gas_particles
    stellar.stop()

    distance = 1 | units.RSun
    sph_secondary.x += distance
    sph_secondary.vx -= 1.7*numpy.sqrt(2*constants.G*stars.mass.sum()/distance)

    sph_particles = Particles()
    sph_particles.add_particles(sph_primary)
    #sph_particles.add_particles(sph_secondary)
    sph_particles.move_to_center()

    converter = nbody_system.nbody_to_si(1|units.hour, 1|units.RSun)

    hydrodynamics = Gadget2(converter)
    hydrodynamics.gas_particles.add_particles(sph_particles)
    hydrodynamics.evolve_model(10.0|units.hour)
    hydrodynamics.gas_particles.copy_values_of_attributes_to(["density", "u",
                                                              "pressure"],
                                                             sph_particles)
    hydrodynamics.stop()

    print "N all=", len(sph_particles)
    clumps = find_clumps(sph_particles, converter)
    #sph_particles = clumps[0]
    print "N blob=", len(sph_particles)
    #plot_clumps(clumps)
    #sph_merger = sph_particles[0]

    print "convert SPH to stellar model"
    merged_star = convert_SPH_to_stellar_model(sph_particles)

    print "initiate stellar evolution model"
    #stellar_evolution = MESA()
    stellar_evolution = EVtwin(redirect="none")    
    stellar_evolution.new_particle_from_model(merged_star, 0.0|units.Myr)
    print "star:", stellar_evolution.particles
    print "evolve star"
    #stellar_evolution.evolve_model(tend)
    while stellar_evolution.model_time<(tend-tcoll):
        stellar_evolution.evolve_model()
        time.append(stellar_evolution.model_time)
        mass.append(stellar_evolution.particles[0].mass)
        radius.append(stellar_evolution.particles[0].radius)
        temperature.append(stellar_evolution.particles[0].temperature)
        luminosity.append(stellar_evolution.particles[0].luminosity)
        print "Time=", time[-1], mass[-1], radius[-1], \
              temperature[-1].in_(units.K), luminosity[-1].in_(units.LSun)

    print stellar_evolution.particles
    stellar_evolution.stop()
    
    return time, mass, radius, temperature, luminosity, n_normal

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", 
                      default = "hydro_triple_gas.hdf5",
                       help="input filename [%default]")
    result.add_option("--tcoll", unit=units.Myr,
                      dest="tcoll", type="float", 
                      default = 50|units.Myr,
                      help="evolution time scale [%default]")
    result.add_option("--tend", unit=units.Myr,
                      dest="tend", type="float", 
                      default = 423|units.Myr,
                      help="evolution time scale [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float", 
                      default = 3|units.MSun,
                      help="stellar mass [%default]")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float", 
                      default = 1|units.MSun,
                      help="stellar mass [%default]")
    return result

if __name__ == "__main__":

    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.RSun,
                                             units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments = new_option_parser().parse_args()
    Mprim = o.Mprim
    Msec = o.Msec
    tend = o.tend
    tcoll = o.tcoll
    
    x_label = "T [K]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True,
                          xsize=14, ysize=10)
    color = get_distinct(4)
    pyplot.xlim(5.e+4, 1.e+3)

    print "Evolve single star of M=", (Mprim).in_(units.MSun)
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(Mprim, tend)

    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[0], lw=2)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[0], s=150, marker="^")

    print "Evolve single star of M=", (Mprim).in_(units.MSun)+(0.2|units.MSun)
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(Mprim+(0.2|units.MSun), tend)

    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[1], lw=2)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[1], s=150, marker="^")

    print "Evolve single star of M=", \
        	(Mprim+Msec).in_(units.MSun) + (0.4|units.MSun)
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(Mprim+(0.4|units.MSun), tend)

    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[3], lw=2)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[3], s=150, marker="^")
    
    tms = 0 |units.Myr
    for i in range(len(stp)):
        if stp[i]>=2 | units.stellar_type:
            tms = time[i]
    if tms <= 1|units.Myr:
        tms = 10|units.Myr
    print "Main-sequence age:", tms.in_(units.Myr)
    tend = tms
    print "Main sequence lifetime of star=", tms.in_(units.Myr)

    #tcoll = 0.5*tms
    time, mass, radius, temperature, luminosity, n \
        = merge_two_stars_sph_and_evolve(o.Mprim, o.Msec, tcoll, o.tend)
    
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[2], s=150, marker="^")
    pyplot.plot(temperature[:n].value_in(units.K),
                luminosity[:n].value_in(units.LSun),
                c=color[2])
    pyplot.scatter(temperature[n-1].value_in(units.K),
                   luminosity[n-1].value_in(units.LSun),
                   c=color[2], s=150, marker="o")
    pyplot.plot(temperature[n:].value_in(units.K),
                luminosity[n:].value_in(units.LSun),
                c=color[2])
    pyplot.scatter(temperature[n+1].value_in(units.K),
                   luminosity[n+1].value_in(units.LSun),
                   c=color[2], s=150, marker="o")

    pyplot.show()

