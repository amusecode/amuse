import numpy
from matplotlib import pyplot

from amuse.lab import *
from amuse.community.mmams.interface import MakeMeAMassiveStarInterface
from amuse.community.mmams.interface import MakeMeAMassiveStar
from amuse.couple.collision_handler import CollisionHandler

from prepare_figure import single_frame
from distinct_colours import get_distinct

default_options = dict()

def evolve_single_star(mass, tend):
    star = Particles(1)
    star.mass = mass
    stellar_evolution = MESA()
    stellar_evolution.particles.add_particles(star)
    time = [] | units.Myr
    stellar_type = [] 
    mass = [] | units.MSun
    radius = [] | units.RSun
    temperature = [] | units.K
    luminosity = [] | units.LSun
    while stellar_evolution.model_time < tend:
        stellar_evolution.evolve_model()
        time.append(stellar_evolution.model_time)
        stellar_type.append(stellar_evolution.particles[0].stellar_type)
        mass.append(stellar_evolution.particles[0].mass)
        radius.append(stellar_evolution.particles[0].radius)
        temperature.append(stellar_evolution.particles[0].temperature)
        luminosity.append(stellar_evolution.particles[0].luminosity)
        print "Time=", time[-1], stellar_type[-1], mass[-1], radius[-1], \
              temperature[-1].in_(units.K), luminosity[-1].in_(units.LSun)
        if stellar_type[-1] >= 4 | units.stellar_type:
            break

    stellar_evolution.stop()
    return time, stellar_type, mass, radius, temperature, luminosity
    
def print_stars(stellar_evolution):
    print "Primary:   Time=", stellar_evolution.model_time.in_(units.Myr), \
        stellar_evolution.particles[0].mass.in_(units.MSun), \
        stellar_evolution.particles[0].radius.in_(units.RSun), \
        stellar_evolution.particles[0].temperature.in_(units.K), \
        stellar_evolution.particles[0].luminosity.in_(units.LSun)
    print "Secondary: Time=", stellar_evolution.model_time.in_(units.Myr), \
        stellar_evolution.particles[1].mass.in_(units.MSun), \
        stellar_evolution.particles[1].radius.in_(units.RSun), \
        stellar_evolution.particles[1].temperature.in_(units.K), \
        stellar_evolution.particles[1].luminosity.in_(units.LSun)

###BOOKLISTSTART1###
def merge_two_stars_and_evolve(Mprim, Msec, tcoll, tend):
    stars = Particles(2)
    stars.mass = [Mprim.value_in(units.MSun),
                  Msec.value_in(units.MSun)] | units.MSun
    stellar_evolution = MESA()
    stellar_evolution.particles.add_particles(stars)
    while stellar_evolution.model_time < tcoll:
        stellar_evolution.evolve_model()
        print_stars(stellar_evolution)
    n_shell = min(stellar_evolution.particles[0].get_number_of_zones(),
                  stellar_evolution.particles[1].get_number_of_zones())
    merger_code = MakeMeAMassiveStar(**default_options)
    merger_code.parameters.target_n_shells = n_shell
    merger_code.parameters.dump_mixed_flag = True
    merger_code.parameters.do_shock_heating_flag = True
    merger_code.commit_parameters()
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
    handler = CollisionHandler(merger_code,
                               stellar_evolution_code = stellar_evolution)
    merger_product = handler.handle_collision(stellar_evolution.particles[0],
                                              stellar_evolution.particles[1])
    merged = stellar_evolution.particles[0]
###BOOKLISTSTOP2###

    print "Stars merged:", merged
    time = [] | units.Myr
    stellar_type = []
    mass = [] | units.MSun
    radius = [] | units.RSun
    temperature = [] | units.K
    luminosity = [] | units.LSun
    
###BOOKLISTSTART3###
    stellar_evolution.evolve_model(keep_synchronous=True)
    p = stellar_evolution.particles[0]
    while stellar_evolution.model_time < tend:
        stellar_evolution.evolve_model()
###BOOKLISTSTOP3###

        time.append(stellar_evolution.model_time)
        stellar_type.append(p.stellar_type)
        mass.append(p.mass)
        radius.append(p.radius)
        temperature.append(p.temperature)
        luminosity.append(p.luminosity)
        
###BOOKLISTSTART4###
        print "Time=", stellar_evolution.model_time, p.stellar_type, \
            p.mass, p.radius, p.temperature.in_(units.K), \
            p.luminosity.in_(units.LSun)
        if p.stellar_type >= 4 | units.stellar_type:
            break
    merger_code.stop()
    stellar_evolution.stop()
###BOOKLISTSTOP4###

    return time, stellar_type, mass, radius, temperature, luminosity

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--tcoll", unit=units.Myr,
                      dest="tcoll", type="float",
                      default = 150|units.Myr,
                      help="moment of collision [%default]")
    result.add_option("--tend", unit=units.Myr,
                      dest="tend", type="float",
                      default = 2|units.Gyr,
                      help="evolution after the collision [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",
                      default = 3|units.MSun,
                      help="Primary ZAMS mass [%default]")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",
                      default = 1|units.MSun,
                      help="Secondary ZAMS mass [%default]")

    return result

if __name__ in ('__main__','__plot__'):

    # High-level structure of merge_two_stars_and_evolve.py and
    # merge_two_stars_sph_evolve.py are designed to be identical.
    
    set_printing_strategy("custom", #nbody_converter = converter, 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    Mprim = o.Mprim
    Msec = o.Msec
    tend = o.tend
    tcoll = o.tcoll

    color = get_distinct(4)  # 0 = cyan, 1 = red, 2 = mustard, 3 = green

    x_label = "T [K]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True,
                          xsize=14, ysize=10)
    pyplot.xlim(2.e+4, 3.e3)
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
    print "Main-sequence lifetime =", tms.in_(units.Myr)
    
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
        = merge_two_stars_and_evolve(Mprim, Msec, tcoll, tend)
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

    save_file = 'merge_two_stars_and_evolve.pdf'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
