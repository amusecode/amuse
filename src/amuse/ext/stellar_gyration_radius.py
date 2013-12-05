"""
   Calculate the radius of gyration for a star
"""
import sys
import numpy
from amuse.lab import *

HeWhiteDwarf = 10 | units.stellar_type
Hertzsprung_gap = 2 | units.stellar_type
First_Asymptotic_Giant_Branch = 5 | units.stellar_type
Second_Asymptotic_Giant_Branch = 6 | units.stellar_type

def calculate_gyration_radius(star):
    I = moment_of_inertia(star)
    k2 = I/(star.mass*star.radius**2)
    return k2

def get_mass_profile(star):
#    if hasattr(star, "get_mass_profile"):
#        mass_profile = star.get_mass_profile()* star.mass
#    else:
    radius_profile = star.get_radius_profile()
    density_profile = star.get_density_profile()
    radii_cubed = radius_profile**3
    radii_cubed.prepend(0|units.m**3)
    mass_profile = (4.0/3.0 * numpy.pi) * density_profile * (radii_cubed[1:] - radii_cubed[:-1])
    print "Derived mass profile from density and radius."
    return mass_profile

def moment_of_inertia(star):
    #Moment of inertia of the Sun: (I/MR^2) = 0.059 
    radius_profile = star.get_radius_profile()
    density_profile = star.get_density_profile()
    I = zero
    dr = radius_profile[1:]-radius_profile[:-1]
    I = density_profile[:-1] * radius_profile[:-1]**4 * dr
    I = I.sum() * (8*numpy.pi/3.) 
    return I

def main(Mstar, z):

#    stellar = MESA()
    stellar = EVtwin()
    stellar.parameters.metallicity = z
    star = stellar.particles.add_particle(Particle(mass=Mstar))

    while star.stellar_type<Second_Asymptotic_Giant_Branch:
        stellar.evolve_model()
        k2 = calculate_gyration_radius(star)
        gamma = 0
        vcrit = (constants.G*star.mass*(1-gamma)/star.radius).sqrt()
        Omega = (vcrit/star.radius)
        J = star.mass *k2 *star.radius**2 * Omega
        print "Star: t=", star.age, "r=", star.radius, "k2=", k2, "J=", J.in_(units.MSun*units.cm**2/units.s), "type=", star.stellar_type

    stellar.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="Mstar", type="float",default = 10.|units.MSun,
                      help="stellar mass [1] %unit")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", 
                          preferred_units = 
                          [units.MSun, units.RSun, units.Myr], 
                          precision = 5, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

