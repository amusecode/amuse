"""
   Visualization for simple N-body integration.
   Reads particle set from file (stellar.hdf5) and prints subsequent frames.
"""
# import sys
# import numpy
from matplotlib import pyplot
from amuse.plot import scatter
from amuse.lab import *
from optparse import OptionParser


def main(filename="stellar.hdf5"):
    Tmax = 100000  # | untis.K
    Tmin = 1000  # | untis.K
    Lmax = 1.e+6  # | untis.LSun
    Lmin = 0.01  # | untis.LSun
    pyplot.ion()
    filename = "nbody.hdf5"
    stars = read_set_from_file(filename, 'hdf5')
    m = 1 + 3.0 * stars.mass / min(stars.mass)
    lim = 2 * max(max(stars.x).value_in(stars.x.unit),
                  stars.center_of_mass().length().value_in(stars.x.unit))
    for si in reversed(list(stars.iter_history())):
        c = si.temperature.value_in(units.K)
        scatter(si.temperature.value_in(units.K),
                si.luminosity.value_in(units.LSun), s=m, c=c,
                cmap=pyplot.cm.hsv)
        pyplot.xlabel("T [K]")
        pyplot.ylabel("L [$L_\odot$]")
        pyplot.xlim(Tmax, Tmin)
        pyplot.ylim(Lmin, Lmax)
        pyplot.loglog()
        pyplot.draw()
        pyplot.cla()
    pyplot.show()


def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default="stellar.hdf5",
                      help="output filename [stellar.hdf5]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)
