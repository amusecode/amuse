import numpy
from matplotlib import pyplot

from amuse.lab import *
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model

from prepare_figure import single_frame
from distinct_colours import get_distinct

def plot_clumps(groups):
    number_of_particles_in_group = []
    fraction_of_mass_in_group =  []
    print "N =", len(groups)
    ci = ['r', 'b', 'g', 'k']
    figure = pyplot.figure(figsize=(6,6))
    i = 0
    alpha = 1
    sizes = 50
    for group in groups:
        pyplot.scatter(group.x.value_in(units.RSun),
                       group.y.value_in(units.RSun),
                       sizes, ci[i], edgecolors = "none", alpha = alpha)
        i += 1
    
    pyplot.xlabel('x (RSun)')
    pyplot.ylabel('y (RSun)')
    size = 100.
    pyplot.xlim(-size, size)
    pyplot.ylim(-size, size)
    pyplot.axis('equal')
    pyplot.savefig('clumps')
    pyplot.show()

def find_clumps(particles, unit_converter):

    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()

    mean_densty = hop.particles.density.mean() 
    hop.parameters.peak_density_threshold = mean_densty
    hop.parameters.saddle_density_threshold = 0.99*mean_densty
    hop.parameters.outer_density_threshold = 0.01*mean_densty

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
    stars.mass = [Mprim.value_in(units.MSun),
                  Msec.value_in(units.MSun)] | units.MSun
    stellar_evolution = MESA()
    stellar_evolution.particles.add_particles(stars)

    # Evolve the stars to tcoll.
    
    while stellar_evolution.model_time < tcoll:
        stellar_evolution.evolve_model()
        print "Time=", stellar_evolution.model_time, \
              stellar_evolution.particles[0].stellar_type, \
              stellar_evolution.particles[0].mass, \
              stellar_evolution.particles[0].radius, \
              stellar_evolution.particles[0].temperature.in_(units.K), \
              stellar_evolution.particles[0].luminosity.in_(units.LSun)

    print stars

    # Convert to SPH particles (nsph per solar mass).

    nsph = 100
    Nprim = int(nsph*stellar_evolution.particles[0].mass.value_in(units.MSun))
    mgas = stellar_evolution.particles[0].mass/Nprim
    Nsec = int(stellar_evolution.particles[1].mass/mgas)
    print "N gas =", Nprim, Nsec
    sph_primary = convert_stellar_model_to_SPH(
        stellar_evolution.particles[0],
        Nprim, 
        seed=12345
    ).gas_particles
    sph_secondary = convert_stellar_model_to_SPH(
        stellar_evolution.particles[1], 
        Nsec, 
        seed=12345
    ).gas_particles
    stellar_evolution.stop()

    # Merge the stars using SPH.
    
    distance = 2 | units.RSun
    sph_secondary.x += distance
    vx = numpy.sqrt(2*constants.G*stars.mass.sum()/distance)
    sph_secondary.vx -= vx
    print 'distance =', distance, 'vx =', vx.in_(units.kms),
    print 'd/v =', (distance/vx).in_(units.hour)
    
    sph_particles = Particles()
    sph_particles.add_particles(sph_primary)
    sph_particles.add_particles(sph_secondary)
    sph_particles.move_to_center()
    
    converter = nbody_system.nbody_to_si(1|units.hour, 1|units.RSun)
    hydrodynamics = Gadget2(converter)
    hydrodynamics.gas_particles.add_particles(sph_particles)
    tf = 10.|units.hour
    hydrodynamics.evolve_model(tf)
    hydrodynamics.gas_particles.copy_values_of_attributes_to(["x", "y", "z",
                                                              "vx", "vy", "vz",
                                                              "density",
                                                              "u", "pressure"],
                                                             sph_particles)
    hydrodynamics.stop()

    # Convert back to a stellar model.
    
    print "N all =", len(sph_particles)
    sph_particles.move_to_center()
    clumps = find_clumps(sph_particles, converter)
    sph_particles = clumps[0]
    print "N blob =", len(sph_particles)
    #plot_clumps(clumps)
    print "convert SPH to stellar model"
    sph_particles.move_to_center()
    merged_star = convert_SPH_to_stellar_model(sph_particles)

    # Evolve the stellar model.
    
    print "initiate stellar evolution model"
    #stellar_evolution = MESA(redirect="none")    
    stellar_evolution = EVtwin(redirect="none")    
    stellar_evolution.new_particle_from_model(merged_star, 0.0|units.Myr)
    print "star:", stellar_evolution.particles
    print "evolve star"
    time = [] | units.Myr
    mass = [] | units.MSun
    radius = [] | units.RSun
    temperature = [] | units.K
    luminosity = [] | units.LSun
    stellar_type = [] 
    while stellar_evolution.model_time < (tend-tcoll):
        try:
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
        except:
            print 'Code crashed at time', stellar_evolution.model_time
            break

    stellar_evolution.stop()
    return time, stellar_type, mass, radius, temperature, luminosity

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--tcoll", unit=units.Myr,
                      dest="tcoll", type="float", 
                      default = 150|units.Myr,
                      help="evolution time scale [%default]")
    result.add_option("--tend", unit=units.Myr,
                      dest="tend", type="float", 
                      default = 2|units.Gyr,
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

    # High-level structure of merge_two_stars_and_evolve.py and
    # merge_two_stars_sph_evolve.py are designed to be identical.
    
    set_printing_strategy("custom", #nbody_converter = converter, 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")
    
    o, arguments = new_option_parser().parse_args()
    Mprim = o.Mprim
    Msec = o.Msec
    tend = o.tend
    tcoll = o.tcoll

    color = get_distinct(4)  # 0 = cyan, 1 = red, 2 = mustard, 3 = green
    
    x_label = "T [K]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True,
                          xsize=14, ysize=10)
    pyplot.xlim(2.e+4, 3.e+3)
    pyplot.ylim(20., 2.e+3)

    print "Evolve single star of mass", Mprim.in_(units.MSun)
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(Mprim, tend)
    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[1], lw=2, zorder=1)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[1], s=150, marker='^',
                   edgecolor='k', zorder=2)

    tms = 0 |units.Myr
    for i in range(len(stp)):
        if stp[i] < 2 | units.stellar_type:
            tms = time[i]
    if tms <= 1|units.Myr:
        tms = 10|units.Myr
    print "Main sequence lifetime =", tms.in_(units.Myr)
    tcoll = 0.5*tms
    icoll = 0
    for i in range(len(stp)):
        if time[i] <= tcoll:
            icoll = i
    pyplot.scatter(temperature[icoll].value_in(units.K),
                   luminosity[icoll].value_in(units.LSun),
                   c=color[2], s=150, marker='o',
                   edgecolor='k', zorder=2)

    print "Evolve single star of mass", (Mprim+Msec).in_(units.MSun)
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(Mprim+Msec, tend)
    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[0], lw=2, zorder=1)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[0], s=150, marker='^',
                   edgecolor='k', zorder=2)

    print "Evolve two single stars and collide at", tcoll.in_(units.Myr)
    time, stp, mass, radius, temperature, luminosity \
        = merge_two_stars_sph_and_evolve(Mprim, Msec, tcoll, tend)
    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[2], ls="--", lw=3, zorder=1)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[2], s=150, marker='o',
                   edgecolor='k', zorder=3)

    Mmerger = mass[0]
    print "Evolve single star of mass", Mmerger
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(Mmerger, tend)
    pyplot.plot(temperature.value_in(units.K),
                luminosity.value_in(units.LSun),
                c=color[3], lw=2, zorder=1)
    pyplot.scatter(temperature[0].value_in(units.K),
                   luminosity[0].value_in(units.LSun),
                   c=color[3], s=150, marker='^',
                   edgecolor='k', zorder=2)

    ax = pyplot.gca()
    ax.tick_params(axis='both', which='both', direction='in')

    save_file = 'merge_two_stars_sph_evolve.pdf'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()

