"""
   Calculate the response of a star as a result of mass loss.
"""
import numpy
from amuse.lab import *

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.RSun, units.Myr], 
                      precision = 6, prefix = "", 
                      separator = " [", suffix = "]")

def calculate_core_mass(star, H_abundance_limit=1.0e-9):
    number_of_zones = star.get_number_of_zones()
    composition = star.get_chemical_abundance_profiles(number_of_zones = number_of_zones)
    index = (composition[0]>H_abundance_limit).nonzero()[0][0] 
    mass = (star.get_cumulative_mass_profile(number_of_zones = number_of_zones) * star.mass)[index]
    return mass[0]

def calculate_zeta(star, z, dmdt= -1*(1.|units.MJupiter)/(1|units.yr)):
    stellar = MESA()
    stellar.parameters.metallicity = z
    stellar.particles.add_particles(star)
    stellar.commit_particles()
    rold = star.radius
    star.mass_change = dmdt
    star.time_step = 1|units.yr
    dm = dmdt*star.time_step
    stellar.particles.evolve_one_step()
    rnew = stellar.particles[0].radius
    dlnr = (rnew-rold)/rold
    dlnm = dm/star.mass
    zeta = dlnr/dlnm
#    if rnew>rold:
#        dlnr = numpy.log((rnew-rold).value_in(units.RSun))
#    else:
#        dlnr = -1*numpy.log((rold-rnew).value_in(units.RSun))
#    zeta = dlnr/numpy.log(-dm.value_in(units.MSun))
#    print "zeta=", (rnew-rold).value_in(units.RSun), -dm.value_in(units.MSun), zeta
    stellar.stop()
    return zeta

def main(Mstar = 1.| units.MSun, dmdt=-0.01| (units.MSun/units.yr), z=0.02):
    stellar = MESA()
    stellar.parameters.metallicity = z
    stellar.particles.add_particles(Particles(mass=Mstar))
    channel_to_framework = stellar.particles.new_channel_to(bodies)
    while stellar.particles[0].stellar_type<Second_Asymptotic_Giant_Branch:
        stellar.particles.evolve_one_step()
        channel_to_framework.copy_attributes(["age", "mass", "radius", "stellar_type"])
        star = stellar.particles.copy_to_memory() 
        zeta = calculate_zeta(star, z, dmdt)
        Mcore = calculate_core_mass(stellar.particles[0])
        print "Zeta=", zeta[0], bodies[0].age, bodies[0].mass, bodies[0].radius, Mcore, dmdt, bodies[0].stellar_type

    stellar.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="Mstar", type="float",default = 1.|units.MSun,
                      help="stellar mass [%default]")
    result.add_option("--dmdt", unit=units.MSun/units.yr,
                      dest="dmdt", type="float",
                      default = -0.01|(units.MSun/units.yr),
                      help="dmdt [%default]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

