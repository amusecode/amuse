"""
   Initialize two stars to a ertain age and merge them using MMAMS
"""
import sys
import numpy
from amuse.lab import *
from amuse.plot import plot, xlabel, ylabel
from matplotlib import pyplot 
    
def merge_two_stars(Mprim, Msec, t_coll):
    bodies = Particles(mass=[Mprim, Msec] |units.MSun)
        
    stellar = MESA()
    primary = stellar.particles.add_particles(bodies[0].as_set())
    secondary = stellar.particles.add_particles(bodies[1].as_set())

    stellar.evolve_model(t_coll)

    print "Pre merger:\n", stellar.particles
    stellar.merge_colliding(primary.copy(), secondary.copy(), MakeMeAMassiveStar,
        dict(), dict(target_n_shells_mixing = 2000), return_merge_products=["se"])
    print "Post merger:\n", stellar.particles

#    composition = stellar.particles[0].get_chemical_abundance_profiles()
    radius = stellar.particles[0].get_radius_profile()
    rho    = stellar.particles[0].get_density_profile()
    print radius
    print rho
    stellar.stop()
    return radius, rho

def plot_density_profile(radius, rho):
    plot(radius.in_(units.RSun), rho)
    pyplot.xlabel("$R$ [$R_\odot$]")
    pyplot.ylabel("density [$g/cm^3$]")
    pyplot.semilogy()
    pyplot.show()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", dest="Mprim", type="float",default = 5,
                      help="Mass of the primary star [%default] MSun")
    result.add_option("-m", dest="Msec", type="float",default = 3,
                      help="Mass of the secondary star [%default] MSun")
    result.add_option("-t", unit=units.Myr, dest="t_coll", type="float", default = 1.0|units.Myr,
                      help="end time of the simulation [%default] Myr")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    radius, rho = merge_two_stars(**o.__dict__)
    plot_density_profile(radius, rho)
