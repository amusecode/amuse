import sys
import numpy
from matplotlib import pyplot 
from amuse.lab import *
from amuse.plot import plot
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model
from prepare_figure import *
from distinct_colours import get_distinct

def return_evolved_star_hydro(mass, time, Nsph):
    star =  Particle(mass=mass)
    stellar = EVtwin()
    star = stellar.particles.add_particle(star)
    stellar.evolve_model(time)
    Nsph = Nsph * int(mass.value_in(units.MSun))
    star_in_sph = convert_stellar_model_to_SPH(star, Nsph).gas_particles
    stellar.stop()
    return star_in_sph

def merge_two_stars_sph(Mprim, Msec, t_coll, Nsph, opening_angle):
    primary_in_sph = return_evolved_star_hydro(Mprim, t_coll,
                                               int(Nsph*Mprim/(1.|units.MSun)))
#    primary_in_sph = relax_sph_realization(primary_in_sph)
    secondary_in_sph = return_evolved_star_hydro(Msec, t_coll,
                                                 int(Nsph*Msec/(1.|units.MSun)))
#    secondary_in_sph = relax_sph_realization(secondary_in_sph)
    R = primary_in_sph.x.max() + secondary_in_sph.x.max()
    M = primary_in_sph.mass.sum() + secondary_in_sph.mass.sum()
    secondary_in_sph.x += 0.8*R
    secondary_in_sph.y += 0.6*R
    secondary_in_sph.vx -= (constants.G*M/R).sqrt()
        
    converter=nbody_system.nbody_to_si(Mprim, 1.0|units.AU)
    hydro = Gadget2(converter, number_of_workers=4)
#    hydro = Fi(converter)
    hydro.parameters.opening_angle = opening_angle

    """
    print "Opening criterion:", hydro.parameters.opening_angle
    print "Opening criterion:", hydro.get_gdgop()
    print "Opening criterion:", hydro.parameters.opening_angle
    print "Opening criterion:", hydro.get_gdgop()
    """
    
    hydro.gas_particles.add_particles(primary_in_sph)
    hydro.gas_particles.add_particles(secondary_in_sph)
    hydro.evolve_model(2.0|units.hour)
    hydro.gas_particles.new_channel_to(primary_in_sph).copy()
    hydro.gas_particles.new_channel_to(secondary_in_sph).copy()
    hydro.stop()
    return primary_in_sph, secondary_in_sph

def relax_sph_realization(sph_star):

    dynamical_timescale = sph_star.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.RSun)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.gas_particles.add_particles(sph_star)

    to_hydro = sph_star.new_channel_to(hydro.gas_particles)
    to_framework = hydro.gas_particles.new_channel_to(sph_star)

    ts_factor = 2.5
    t_end = ts_factor * sph_star.dynamical_timescale(mass_fraction=0.9)
    n_steps = ts_factor * 100
    velocity_damp_factor = 1.0 - (ts_factor*2*numpy.pi)/n_steps
    dt = t_end/float(n_steps)
    time = 0|units.day
    while time < t_end:
        time += dt
        hydro.evolve_model(time)
        hydro.gas_particles.velocity = velocity_damp_factor * hydro.gas_particles.velocity
    to_framework.copy()
    hydro.stop()
    return sph_star

def merge_and_plot_distribution(Mprim, Msec, t_coll, Nsph, opening_angle, c, label):
    p, s = merge_two_stars_sph(Mprim, Msec, t_coll, Nsph, opening_angle)
    merger = ParticlesSuperset([p, s])
    com = merger.center_of_mass()
    merger.r = ((merger.x-com[0])**2 + (merger.y-com[1])**2 + (merger.z-com[2])**2).sqrt()
    merger = merger.sorted_by_attributes("r")
    n = []
    m = 0
    mi =  1.0*(Mprim+Msec).value_in(units.MSun)/len(merger)
    for i in range(len(merger.r)):
        m += mi
        n.append(m)
    pyplot.plot(merger.r.value_in(units.RSun), n, c=c, label=label)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",default = 10|units.MSun,
                      help="Mass of the primary star [%default] MSun")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",default = 1|units.MSun,
                      help="Mass of the secondary star [%default] MSun")
    result.add_option("-N", 
                      dest="Nsph", type="int",default = 100,
                      help="Number of sph particles per MSun [%default]")
    result.add_option("-t", unit=units.Myr, 
                      dest="t_coll", type="float", default = 0.01|units.Myr,
                      help="end time of the simulation [%default] Myr")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    colors = get_distinct(4)
    figure = pyplot.figure(figsize=(16, 12))
    pyplot.xlabel("R [R$_\odot$]")
    pyplot.ylabel("$M_{<r}$")
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    merge_and_plot_distribution(o.Mprim, o.Msec, o.t_coll, 10, 0.5, colors[0], "$N_{sph} = 10/M_\odot$")
    merge_and_plot_distribution(o.Mprim, o.Msec, o.t_coll, 100, 0.5, colors[1], "$N_{sph} = 100/M_\odot$")
    merge_and_plot_distribution(o.Mprim, o.Msec, o.t_coll, 1000, 0.5, colors[2], "$N_{sph} = 10^3/M_\odot$")
    merge_and_plot_distribution(o.Mprim, o.Msec, o.t_coll, 10000, 0.5, colors[3], "$N_{sph} = 10^4/M_\odot$")
#    pyplot.show()
#    pyplot.semilogx()
    pyplot.xlim(0, 20)
    pyplot.ylim(0, 12)
    pyplot.legend(loc=4, fontsize=24)
    pyplot.savefig("stellar_merger_convergence")
