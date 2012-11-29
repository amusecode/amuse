"""
   Simple visualization for N-body integration.
   Reads particle set from file (nbody.hdf5) and prints frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from amuse.units.optparse import OptionParser
from time import sleep

def main(filename, lim=-1|units.parsec):
    pyplot.ion()
    mc = read_set_from_file(filename, "hdf5") # particles
    if lim<=zero:
        lim = max(mc.x).value_in(lim.unit)
#    lim = max(stars.x).value_in(stars.x.unit)
    time = 0
    for si in mc.history:
#        time = si.get_timestamp()
        print si
        m = 1
        pyplot.title("Cluster at t="+str(time))
#        print "time = ", time, mc[0].position.value_in(units.parsec)
#        scatter(si.x, si.y, s=m)
        scatter(si.x.as_quantity_in(lim.unit), si.y.as_quantity_in(lim.unit), s=m)
#        scatter(si.x, si.y, s=m)
        #scatter(mc[0].x, mc[0].y, c='r', s=100.)
        xlabel("X")
        ylabel("Y")
        if lim>zero:
            pyplot.xlim(-lim.value_in(lim.unit), lim.value_in(lim.unit))
            pyplot.ylim(-lim.value_in(lim.unit), lim.value_in(lim.unit))
        pyplot.draw()
#        sleep(1)
        pyplot.cla()
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity.hdf5",
                      help="output filename [gravty.hdf5]")
    result.add_option("-l", unit=units.parsec, 
                      dest="lim", type="float", default = 1000|units.parsec,
                      help="axis length [1000] %unit")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


