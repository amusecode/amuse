import sys
import numpy
from matplotlib import pyplot 
from amuse.lab import *
from amuse.plot import plot
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model

###BOOKLISTSTART1###
def return_evolved_star_hydro(mass, time, Nsph):
    star =  Particle(mass=mass)
    stellar = EVtwin()
    star = stellar.particles.add_particle(star)
    stellar.evolve_model(time)
    Nsph = Nsph * int(mass.value_in(units.MSun))
    star_in_sph = convert_stellar_model_to_SPH(star, Nsph).gas_particles
    stellar.stop()
    return star_in_sph
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def merge_two_stars_sph(Mprim, Msec, t_coll, Nsph):
    primary_in_sph = return_evolved_star_hydro(Mprim, t_coll, Nsph)
    primary_in_sph = relax_sph_realization(primary_in_sph)
    secondary_in_sph = return_evolved_star_hydro(Msec, t_coll, Nsph)
    secondary_in_sph = relax_sph_realization(secondary_in_sph)
    R = primary_in_sph.x.max() + secondary_in_sph.x.max()
    M = primary_in_sph.mass.sum() + secondary_in_sph.mass.sum()
    secondary_in_sph.x += 0.8*R
    secondary_in_sph.y += 0.6*R
    secondary_in_sph.vx -= (constants.G*M/R).sqrt()
        
    converter=nbody_system.nbody_to_si(Mprim, 1.0|units.AU)
    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(primary_in_sph)
    hydro.gas_particles.add_particles(secondary_in_sph)
    hydro.evolve_model(2.0|units.hour)
    hydro.gas_particles.new_channel_to(primary_in_sph).copy()
    hydro.gas_particles.new_channel_to(secondary_in_sph).copy()
    hydro.stop()
    return primary_in_sph, secondary_in_sph
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
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
###BOOKLISTSTOP3###

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
    p, s = merge_two_stars_sph(**o.__dict__)
    pyplot.scatter(p.x.value_in(units.RSun), p.y.value_in(units.RSun), c='b')
    pyplot.scatter(s.x.value_in(units.RSun), s.y.value_in(units.RSun), c='r')
    pyplot.show()
