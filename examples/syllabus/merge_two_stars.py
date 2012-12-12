"""
   Initialize two stars to a ertain age and merge them using MMAMS
"""
import sys
import numpy
from amuse.lab import *
    
def merge_two_stars(Mprim, Msec, t_coll):
    bodies = Particles(mass=[Mprim, Msec] |units.MSun)
        
    stellar = MESA()
    primary = stellar.particles.add_particles(bodies[0].as_set())
    secondary = stellar.particles.add_particles(bodies[1].as_set())

    stellar.evolve_model(t_coll)

    print "Pre merger:\n", stellar.particles
    stellar.merge_colliding(primary, secondary, MakeMeAMassiveStar,
                            dict(), dict())
    print "Post merger:\n", stellar.particles

    composition = stellar.particles[0].get_composition_profile()
    print "composition:\n", composition


    stellar.stop()

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
    merge_two_stars(**o.__dict__)
