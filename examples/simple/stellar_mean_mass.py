"""
Calculate the mean mass of a stellar popultion as a function of time.
This routine was used to measure the <m> of a star cluster in order to
incorporate <m>(t) and d<m>/dt in a parametrized cluster evolution code.
"""

import sys
import numpy
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel

from amuse.community.sse.interface import SSE
from amuse.lab import *

from optparse import OptionParser
    
if __name__ in ('__main__', '__plot__'):        

    result = OptionParser()
    result.add_option("-t", dest="t_end", type="float",default = 12000,
                      help="End time for calculation in Myr. [12000 Myr]")
    result.add_option("-d", dest="dt", type="float",default = 10,
                      help="output time step in Myr. [10 Myr]")
    result.add_option("-N", dest="N", type="int",default = 10,
                      help="number of stars to calculate <m>. [10]")
    result.add_option("-m", dest="Mmin", type="float",default = 0.15,
                      help="Minimal mass of the IMF in MSun. [0.15MSun]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="Maximal mass of the IMF in MSun. [100MSun]")
    result.add_option("-x", dest="x_imf", type="float",default = -2.35,
                      help="Slope of the IMF. [-2.35]")
    result.add_option("-v", dest="verbose", action="store_true",default=False,
                      help="verbose output [True]")

    o, arguments  = result.parse_args()
    t_end = o.t_end | units.Myr
    dt = o.dt | units.Myr

    stellar_evolution = SSE()
    stellar_evolution.commit_parameters()

    stars=Particles(o.N)
    Mmin = o.Mmin | units.MSun
    Mmax = o.Mmax | units.MSun
    if o.verbose:
        print "#Selected parameters: "
        print "#\tN=", o.N
        print "#\tIMF=", o.Mmin, "MSun", o.Mmax, "MSun", o.x_imf
        print "#\t t [Myr] \t <m> [MSun] \t\t d<m>/dt [MSun/Myr]"
        
    stars.mass = new_salpeter_mass_distribution(
        o.N, mass_min=Mmin, mass_max=Mmax, alpha=o.x_imf)

    stars = stellar_evolution.particles.add_particles(stars)
    stellar_evolution.commit_particles()
    t = 0 | units.Myr
    mm = stars.mass.sum()/len(stars)
    while t<t_end:
        mm_last = mm
        t += dt
        stellar_evolution.evolve_model(t)
        mm = stars.mass.sum()/len(stars)
        dmm_dt = (mm_last-mm)/dt
        if o.verbose:
            print "t = ", t, "<m>=", mm.as_quantity_in(units.MSun), " d<m>/dt = ", dmm_dt.as_quantity_in(units.MSun/units.Myr)
        else:
            print "\t", t, "\t", mm.as_quantity_in(units.MSun), " \t ", dmm_dt.as_quantity_in(units.MSun/units.Myr)

