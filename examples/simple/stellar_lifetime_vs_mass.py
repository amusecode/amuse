""" 
Calculates the stellar lifetime in a range of masses between 
Mmax and Mmin using SSE (or another stellar evolution code)
and an analytic expression.
"""

import numpy
from optparse import OptionParser
from amuse.lab import *
from amuse.plot import plot
from matplotlib import pyplot as plt
from amuse.community.sse.interface import SSE

se = None

def stellar_remnant_state(star):
    return 10 <= star.stellar_type.value_in(units.stellar_type) and \
        star.stellar_type.value_in(units.stellar_type) < 16

def stellar_lifetime(mZAMS, z=0.02):
    global se
    if se is None:
        se = SSE()
        se.parameters.metallicity = z
    se.particles.add_particle(Particle(mass=mZAMS))
    while not stellar_remnant_state(se.particles[0]):
        se.evolve_model()
    t_end = se.particles[0].age
    tpe = se.particles[0].stellar_type
    se.particles.remove_particle(se.particles[0])
    return t_end

def power_law_fit_to_main_sequence_lifetime(mZAMS):
    return 2 + 1.0E+4/pow(mZAMS.value_in(units.MSun), 2.5) | units.Myr

def main(n=10, mmin=1.0, mmax=100, z=0.02):
    dm = (mmax-mmin)/n
    mZAMS = numpy.arange(mmin, mmax, dm) | units.MSun
    mmin=mmin|units.MSun
    mmax=mmax|units.MSun
    print mZAMS
    t_sse = [] | units.Myr
    t_analytic = [] | units.Myr
    for mi in mZAMS:
        t_sse.append(stellar_lifetime(mi, z))
        t_analytic.append(power_law_fit_to_main_sequence_lifetime(mi))
    plot(mZAMS, t_sse, label="sse")
    plot(mZAMS, t_analytic,label="analytic")
    plt.loglog()
    plt.legend()
    plt.title("comparison between SSE and analytic with z="+str(z))
    plt.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n", type="int",default = 10,
                      help="number of stars")
    result.add_option("-m", dest="mmin", type="float", default = 1.0,
                      help="Minimal mass [1.0] MSun")
    result.add_option("-M", dest="mmax", type="float", default = 100.0,
                      help="Maximal mass [100] MSun")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

