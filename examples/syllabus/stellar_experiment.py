"""
   Calculate the response of a star as a result of mass loss.
"""
from amuse.lab import *

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.RSun, units.Myr], 
                      precision = 6, prefix = "", 
                      separator = " [", suffix = "]")

def calculate_zeta(star, z, dmdt):
    stellar = MESA()
    stellar.parameters.metallicity = z
    stellar.particles.add_particles(star)
    stellar.commit_particles()
    rold = star.radius
    star.mass_change = dmdt
    dm = 0.01*star.mass
    star.time_step = dm/dmdt
    stellar.particles.evolve_one_step()
    rnew = stellar.particles[0].radius
    dlnr = (rnew-rold)/rold
    dlnm = (stellar.particles[0].mass-star.mass)/star.mass
    zeta = dlnr/dlnm

    stellar.stop()
    return zeta

def main(Mstar, z, dmdt):
    stellar = MESA()
    stellar.parameters.metallicity = z
    bodies = Particles(mass=Mstar)
    stellar.particles.add_particles(bodies)
    channel_to_framework = stellar.particles.new_channel_to(bodies)
    copy_argument = ["age", "mass", "radius", "stellar_type"]
    while stellar.particles[0].stellar_type<Second_Asymptotic_Giant_Branch:
        stellar.particles.evolve_one_step()
        channel_to_framework.copy_attributes(copy_argument)
        star = stellar.particles.copy() 
        zeta = calculate_zeta(star, z, dmdt)
        print "Zeta=", zeta[0], bodies[0].age, bodies[0].mass, bodies[0].radius, dmdt, bodies[0].stellar_type

    stellar.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun, dest="Mstar", type="float",default = 1.|units.MSun, help="stellar mass [%default]")
    result.add_option("--dmdt", unit=units.MSun/units.yr, dest="dmdt", type="float", default = -0.01|(units.MSun/units.yr), help="dmdt [%default]")
    result.add_option("-z", dest="z", type="float", default = 0.02, help="metalicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

