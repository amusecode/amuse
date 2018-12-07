import os.path
import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.evrard_test import uniform_unit_sphere
import time

from initialize_sstars import *
from amuse.community.adaptb.interface import Adaptb
from gravity_class import Gravity

#from prepare_figure import single_frame, figure_frame, set_tickmarks
#from distinct_colours import get_distinct

def _plot_orbits(x, y, z):
    from matplotlib import pyplot, rc
    fig = pyplot.figure(figsize=(14,14))
    font = {'size' : 30}
    rc('font', **font)
    pyplot.figaspect(1.0)
    pyplot.scatter(x[0].value_in(units.AU), y[0].value_in(units.AU), 200, lw=0)
    pyplot.plot(x.value_in(units.AU), y.value_in(units.AU), lw=2)
    pyplot.xlabel("$X [AU]$")
    pyplot.ylabel("$Y [AU]$")
    pyplot.xlim(-15000, 15000)
    pyplot.ylim(-30000, 10000)
#    pyplot.show()
    pyplot.savefig("SStars_1Jan2001_orbits")

def plot_orbits(x, y, z):
    from matplotlib import pyplot, rc
    fig = pyplot.figure(figsize=(14,14))
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)
    font = {'size' : 30}
    rc('font', **font)
    pyplot.figaspect(1.0)
    pyplot.scatter(x[0].value_in(units.AU), y[0].value_in(units.AU), 200, lw=0)
    pyplot.plot(x.value_in(units.AU), y.value_in(units.AU), lw=2)
    pyplot.xlabel("$X [AU]$")
    pyplot.ylabel("$Y [AU]$")
    pyplot.xlim(-10000, 10000)
    pyplot.ylim(-10000, 10000)
#    pyplot.show()
    pyplot.savefig("SStars_1Jan2001_orbits")

def main(t_end=1, n_steps=1, filename=None):

    black_hole, stars = initialize_sstars(2001|units.yr, S_name, S_a_arcsec, S_ecc, S_inc, S_omra, S_Omega, S_tperi, S_Period)
    print "N=", len(stars)

    gravity = Gravity(Mercury, [black_hole, stars])
    model_time = gravity.model_time

    model_time = 0 | units.Myr
    dt_diag = t_end/float(n_steps)
    t_diag = model_time
    dt = dt_diag
    x = [] | units.AU
    y = [] | units.AU
    z = [] | units.AU
    while model_time < t_end:
        model_time += dt
        gravity.evolve_model(model_time)
        x.append(gravity.particles.x)
        y.append(gravity.particles.y)
        z.append(gravity.particles.z)

    gravity.stop()
    plot_orbits(x, y, z)    

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 100,
                      help="number of diagnostics time steps [10]")
    result.add_option("-f", dest="filename", default = None,
                      help="write output filename")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 500|units.yr,
#                      dest="t_end", type="float", default = 1678|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
