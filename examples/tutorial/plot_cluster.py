"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser

def main(filename = "nbody.hdf5", lim=3):
    pyplot.ion()
    storage = store.StoreHDF(filename,"r")
    stars = storage.load()
#    lim = 4*stars.center_of_mass().length().value_in(stars.x.unit)
    lim = 2*max(max(stars.x).value_in(stars.x.unit),  stars.center_of_mass().length().value_in(stars.x.unit))
    m =  1 + 3.0*stars.mass/min(stars.mass)
    for si in stars.history:
        time = si.get_timestamp()
        pyplot.title("Cluster at t="+str(time))
        print "time = ", time
        scatter(si.x, si.y, s=m)
        xlabel("X")
        ylabel("Y")
        pyplot.xlim(-lim, lim)
        pyplot.ylim(-lim, lim)
        pyplot.draw()
        pyplot.cla()
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


