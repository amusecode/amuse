"""
   Simple visualization for N-body integration.
   Reads particle set from file (nbody.hdf5) and prints frames.
"""
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store

def main(filename, lim):
    pyplot.ion()
    particles = read_set_from_file(filename, "hdf5") 
    if lim<=zero:
        lim = max(particles.x).value_in(lim.unit)
    time = 0
    for si in particles.history:
        pyplot.title("Cluster at t="+str(time))
        scatter(si.x.as_quantity_in(lim.unit), si.y.as_quantity_in(lim.unit))
        xlabel("X")
        ylabel("Y")
        if lim>zero:
            pyplot.xlim(-lim.value_in(lim.unit), lim.value_in(lim.unit))
            pyplot.ylim(-lim.value_in(lim.unit), lim.value_in(lim.unit))
        pyplot.draw()
        pyplot.cla()
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity.hdf5", help="output filename [gravty.hdf5]")
    result.add_option("-l", unit=units.kpc, 
                      dest="lim", type="float", default = -1|units.kpc, help="axis length [1000] %unit")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


