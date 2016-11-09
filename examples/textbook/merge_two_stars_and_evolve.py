import numpy
from amuse.lab import *
from amuse.community.mmams.interface import MakeMeAMassiveStarInterface, MakeMeAMassiveStar
from amuse.couple.collision_handler import CollisionHandler

from matplotlib import pyplot
from prepare_figure import single_frame
from distinct_colours import get_distinct

default_options = dict()

def get_density_profile(code=MESA, M=1.0|units.MSun, z=0.02, t=2|units.Myr):
    stellar = code()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=M))
    print "Nzones=", stellar.particles.get_number_of_zones()
    stellar.evolve_model(t)
    radius = stellar.particles[0].get_radius_profile()
    rho    = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho

def print_stars(stellar_evolution):
    print "Time= Primary", stellar_evolution.model_time.in_(units.Myr), stellar_evolution.particles[0].mass.in_(units.MSun), stellar_evolution.particles[0].radius.in_(units.RSun), stellar_evolution.particles[0].temperature.in_(units.K), stellar_evolution.particles[0].luminosity.in_(units.LSun)
    print "Time= Secondary", stellar_evolution.model_time.in_(units.Myr), stellar_evolution.particles[1].mass.in_(units.MSun), stellar_evolution.particles[1].radius.in_(units.RSun), stellar_evolution.particles[1].temperature.in_(units.K), stellar_evolution.particles[1].luminosity.in_(units.LSun)

def merge_two_stars(Mprim, Msec, n_shell, tcoll, tend):
        stars = Particles(2)
        stars.mass = [Mprim.value_in(units.MSun), Msec.value_in(units.MSun)] | units.MSun
        
        stellar_evolution = MESA()
        stellar_evolution.particles.add_particles(stars)
        while stellar_evolution.model_time<tcoll:
            stellar_evolution.evolve_model()
            print_stars(stellar_evolution)

        n_shell = min(stellar_evolution.particles[0].get_number_of_zones(), stellar_evolution.particles[1].get_number_of_zones())
        print "n_shells=", n_shell

        instance = MakeMeAMassiveStar(**default_options)
        instance.parameters.target_n_shells = n_shell
        instance.parameters.dump_mixed_flag = True
        instance.parameters.do_shock_heating_flag = True
        instance.commit_parameters()
        instance.particles.add_particles(stars)

        handler = CollisionHandler(instance, stellar_evolution_code = stellar_evolution)
        merger_product = handler.handle_collision(stellar_evolution.particles[0], stellar_evolution.particles[1])
        merged = stellar_evolution.particles[0]
        stellar_evolution.evolve_model(keep_synchronous = True)

        time = [] | units.Myr
        mass = [] | units.MSun
        radius = [] | units.RSun
        temperature = [] | units.K
        luminosity = [] | units.LSun
        while stellar_evolution.model_time<tend:
#        while stellar_evolution.particles[0].radius<24|units.RSun:
            stellar_evolution.evolve_model()
            time.append(stellar_evolution.model_time)
            mass.append(stellar_evolution.particles[0].mass)
            radius.append(stellar_evolution.particles[0].radius)
            temperature.append(stellar_evolution.particles[0].temperature)
            luminosity.append(stellar_evolution.particles[0].luminosity)
            print "Time=", time[-1], mass[-1], radius[-1], temperature[-1].in_(units.K), luminosity[-1].in_(units.LSun)

        rho_profile = merged.get_density_profile()
        radius_profile = merged.get_radius_profile()

        instance.stop()
        stellar_evolution.stop()

        return time, mass, radius, temperature, luminosity

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--tcoll", unit=units.Myr,
                      dest="tcoll", type="float",default = 1|units.Myr,
                      help="moment of collision [%default]")
    result.add_option("--tend", unit=units.Myr,
                      dest="tend", type="float",default = 1|units.Myr,
                      help="evolution after the collision [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",default = 1|units.MSun,
                      help="Primary ZAMS mass [%default]")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",default = 1|units.MSun,
                      help="Secondary ZAMS mass [%default]")
    result.add_option("-n", 
                      dest="n_shell", type="int",default = 1000,
                      help="number of mass shells [%default]")
    return result

def plot_post_collision_star(time, mass, radius, temperature, luminosity):
    pyplot.subplot(2,2,1)
    pyplot.plot(time.value_in(units.Myr), mass.value_in(units.MSun))
    pyplot.subplot(2,2,2)
    pyplot.plot(time.value_in(units.Myr), radius.value_in(units.RSun))
    pyplot.subplot(2,2,3)
    pyplot.plot(temperature.value_in(units.K), luminosity.value_in(units.LSun))
    pyplot.loglog()
    pyplot.show()

if __name__ in ('__main__','__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          precision = 15, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()

    x_label = "T [K]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=True, logy=True, xsize=14, ysize=10)
    color = get_distinct(3)

    Mprim = 1|units.MSun
    Msec = 1|units.MSun
    tend = 100|units.Myr
    n_shell = 1000

    tcoll = 0|units.Myr
    time, mass, radius, temperature, luminosity = Evolve_single_star(Mprim+Msec, tend)
    pyplot.plot(temperature.value_in(units.K), luminosity.value_in(units.LSun), c=color[0])

    tcoll = 50|units.Myr
    time, mass, radius, temperature, luminosity = merge_two_stars(Mprim, Msec, n_shell, tcoll, tend)
    pyplot.plot(temperature.value_in(units.K), luminosity.value_in(units.LSun), c=color[1])
    
    tcoll = tend
    time, mass, radius, temperature, luminosity = merge_two_stars(Mprim, Msec, n_shell, tcoll, tend)
    pyplot.plot(temperature.value_in(units.K), luminosity.value_in(units.LSun), c=color[2])
    
    
#    plot_post_collision_star(time, mass, radius, temperature, luminosity)

    """
    from matplotlib import pyplot
    pyplot.plot(r.value_in(units.RSun), m.value_in(units.g/units.cm**3), c='b')

    r, m = get_density_profile(code=MESA, M=o.Mprim, z=0.02, t=o.tend)
    pyplot.plot(r.value_in(units.RSun), m.value_in(units.g/units.cm**3), ls='--', c='r')

    r, m = get_density_profile(code=MESA, M=o.Msec, z=0.02, t=o.tend)
    pyplot.plot(r.value_in(units.RSun), m.value_in(units.g/units.cm**3), ls='-.', c='b')

    r, m = get_density_profile(code=MESA, M=o.Mprim+o.Msec, z=0.02, t=o.tend)
    pyplot.plot(r.value_in(units.RSun), m.value_in(units.g/units.cm**3), c='g')
    pyplot.show()
    """
