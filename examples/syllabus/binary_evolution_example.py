"""
Generates a grid of binaries with different, primary mass, mass ratio
and separation and evolves these over time.
"""

from amuse.units import units
from amuse.units import quantities
from amuse import datamodel
from amuse.community.seba.interface import SeBa
from amuse.community.bse.interface import BSE
from matplotlib import pyplot

import numpy
import time

def create_double_star(Mprim, Msec, a, e):
    primary_stars   = datamodel.Particles(mass=Mprim)
    secondary_stars = datamodel.Particles(mass=Msec)
        
    stars = datamodel.Particles()
    primary_stars   = stars.add_particles(primary_stars)
    secondary_stars = stars.add_particles(secondary_stars)
        
    double_star = datamodel.Particles(
        semi_major_axis = a,
        eccentricity    = e
        )
    double_star.child1 = list(primary_stars)
    double_star.child2 = list(secondary_stars)
    return double_star, stars

def evolve_double_star(Mprim, Msec, a, e, end_time, n_steps):
    q = Msec/Mprim
    double_star, stars = create_double_star(Mprim, Msec, a, e)
    time = 0|units.Myr
    time_step = end_time/n_steps

#    code = BSE()
    code = SeBa()
    code.particles.add_particles(stars)
    code.binaries.add_particles(double_star)
    
    channel_from_code_to_model_for_binaries = code.binaries.new_channel_to(double_star)
    channel_from_code_to_model_for_stars = code.particles.new_channel_to(stars)
    
    #we evolve in steps of timestep, just to get some feedback
    print "start evolving..."
    t = []
    a = []
    e = []
    while time < end_time:
        time += time_step
        code.evolve_model(time)
#        print "evolved to time: ", time.as_quantity_in(units.Myr), double_star[0].semi_major_axis.as_quantity_in(units.RSun), double_star[0].eccentricity, stars[0].mass, stars[1].mass

        channel_from_code_to_model_for_stars.copy()
        channel_from_code_to_model_for_binaries.copy()
        t.append(time.value_in(units.Myr))
        a.append(double_star[0].semi_major_axis.value_in(units.RSun))
        e.append(double_star[0].eccentricity)
    code.stop()

    pyplot.figure(figsize = (8,8))
    pyplot.title('Binary evolution', fontsize=12)
    pyplot.plot(t, a)
    pyplot.xlabel('time [Myr]')
    pyplot.ylabel('semi major axis (AU)')
#    pyplot.ylabel('eccentricity')
    pyplot.show()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",default = 12|units.MSun,
                      help="primar mass [%defailt]")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",default =  10|units.MSun,
                      help="secondary mass [%defailt]")
    result.add_option("-T", unit=units.Myr,
                      dest="end_time", type="float", default = 100.0 |units.Myr,
                      help="end time of the simulation [%defailt]")
    result.add_option("-a", unit=units.RSun,
                      dest="a", type="float",default =  205|units.RSun,
                      help="orbital separation [%defailt]")
    result.add_option("-e", dest="e", type="float", default = 0,
                      help="orbital eccentricity [%defailt]")
    result.add_option("-n", dest="n_steps", type="float", default = 10,
                      help="number of output steps [%defailt]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    evolve_double_star(**o.__dict__)

