"""
   Calculate the response of a star as a result of mass loss.
"""
from amuse.lab import *

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
MASS_UNIT = units.MSun
LENGTH_UNIT = units.RSun
TIME_UNIT = units.Myr
MASSLOSS_UNIT = units.MSun / units.yr
set_printing_strategy("custom", 
                      preferred_units = [MASS_UNIT, LENGTH_UNIT, TIME_UNIT, MASSLOSS_UNIT], 
                      precision = 6, prefix = "", 
                      separator = " [", suffix = "]")

###BOOKLISTSTART1###
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
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def main(Mstar, z, dmdt):
    stellar = MESA()
    stellar.parameters.metallicity = z
    bodies = Particles(mass=Mstar)
    stellar.particles.add_particles(bodies)
    channel_to_framework = stellar.particles.new_channel_to(bodies)
    copy_argument = ["age", "mass", "radius", "stellar_type"]
    while stellar.particles[0].stellar_type < Second_Asymptotic_Giant_Branch:
        stellar.particles.evolve_one_step()
        channel_to_framework.copy_attributes(copy_argument)
        star = stellar.particles.copy() 
        zeta = calculate_zeta(star, z, dmdt)
        print("Zeta=", zeta[0], bodies[0].age, bodies[0].mass, \
              bodies[0].radius, dmdt, bodies[0].stellar_type)
    stellar.stop()
###BOOKLISTSTOP2###
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=MASS_UNIT, dest="Mstar", type="float",
                      default=1. | MASS_UNIT, help="stellar mass [%default]")
    result.add_option("--dmdt", unit=MASSLOSS_UNIT, dest="dmdt",
                      type="float", default=-0.01 | (MASSLOSS_UNIT),
                      help="dmdt [%default]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metallicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

