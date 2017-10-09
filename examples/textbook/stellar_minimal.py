"""
   Minimal routine for running a stellar evolution code
"""

###BOOKLISTSTART###
from amuse.lab import *

def main(m, z, model_time):
    stellar = MESA()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=m))

    initial_luminosity = stellar.particles.luminosity
    dt = 1 | units.Myr
    while stellar.model_time < model_time:
        stellar.evolve_model(stellar.model_time+dt)

    print "at T=", stellar.model_time.in_(units.Myr), \
        "L(t=0)=", initial_luminosity, \
        ", L (t=", stellar.particles.age.in_(units.Myr), \
        ")=", stellar.particles.luminosity.in_(units.LSun), \
        ", m=", stellar.particles.mass.in_(units.MSun), \
        ", R=", stellar.particles.radius.in_(units.RSun)
        
    stellar.stop()
###BOOKLISTSTOP###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-m", unit= units.MSun,
                      dest="m", type="float", default = 1.0 | units.MSun,
                      help="stellar mass [%default]")
    result.add_option("-t", unit = units.Myr,
                      dest="model_time", type="float", 
                      default = 4700.0|units.Myr,
                      help="end time of the simulation [%default]")
    result.add_option("-z", dest="z", type="float", 
                      default = 0.02, help="metalicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
