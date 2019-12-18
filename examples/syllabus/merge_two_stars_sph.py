"""
   Initialize two stars to a ertain age and merge them using MMAMS
"""
import sys
import numpy
from amuse.lab import *
from amuse.plot import plot, xlabel, ylabel
from matplotlib import pyplot 

from orbital_elements_to_Cartesian import Orbital_elements_to_pos_vel
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model

def convert_star_to_hydro_model(M, t_end):

    return sph_particles

def plot_density_profile(radius, rho):
    plot(radius.in_(units.RSun), rho)
    pyplot.xlabel("$R$ [$R_\odot$]")
    pyplot.ylabel("density [$g/cm^3$]")
#    pyplot.semilogy()
#    pyplot.show()
    
def merge_two_stars_sph(Mprim, Msec, t_coll):

    star =  Particle(mass=Mprim)
    stellar_evolution = EVtwin()
    EVTwin_star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.evolve_model(t_coll)
    print("star=", EVTwin_star)
    EVTwin_radius = EVTwin_star.get_radius_profile()
    EVTwin_rho    = EVTwin_star.get_density_profile()
    
    N_sph = 100*Mprim.value_in(units.MSun)
    primary_star = convert_stellar_model_to_SPH(EVTwin_star, N_sph).gas_particles
    stellar_evolution.stop()

    converter=nbody_system.nbody_to_si(Mprim, 1.0|units.AU)
    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(primary_star)
    channel =  hydro.gas_particles.new_channel_to(primary_star)
    hydro.evolve_model(1.0|units.s)

    channel.copy()
    hydro.stop()

    pyplot.scatter(primary_star.x.value_in(units.AU), primary_star.y.value_in(units.AU))
    pyplot.show()

    new_stellar_model = convert_SPH_to_stellar_model(primary_star)
    stellar_evolution = MESA(redirection="none")
    stellar_evolution.commit_parameters()
    
    stellar_evolution.new_particle_from_model(new_stellar_model, t_coll)
    MESA_star = stellar_evolution.particles[0]
    print("star=", MESA_star)
    MESA_radius = MESA_star.get_radius_profile()
    MESA_rho    = MESA_star.get_density_profile()
    stellar_evolution.stop()
    plot_density_profile(MESA_radius, MESA_rho)
    plot_density_profile(EVTwin_radius, EVTwin_rho)
    pyplot.semilogy()
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",default = 10|units.MSun,
                      help="Mass of the primary star [%default] MSun")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",default = 1|units.MSun,
                      help="Mass of the secondary star [%default] MSun")
    result.add_option("-t", unit=units.Myr, 
                      dest="t_coll", type="float", default = 0.01|units.Myr,
                      help="end time of the simulation [%default] Myr")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    merge_two_stars_sph(**o.__dict__)
