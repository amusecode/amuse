"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
#from optparse import OptionParser
from time import sleep
from amuse.units.optparse import OptionParser
from amuse.plot import sph_particles_plot

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import numpy as np

from prepare_figure import single_frame
#from distinct_colours import get_distinct


def main(filename = "hydro.hdf5", lim=None, image_id=-1, nice_plot=1):
    x_label = 'x'
    y_label = 'y'
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=10)
    stars = read_set_from_file(filename, "hdf5")
    print(stars)
    snapshot_id = 0
    isgas = True
    snapshot_id = 0
    sinks = Particles(0)    
    for si in stars.history:
        if isgas:
            gas = si
            isgas = False
        else:
            sinks = si
            isgas = True
        if not isgas:
            snapshot_id += 1
            time = gas.get_timestamp()
            if image_id<0 or image_id == snapshot_id:
                pyplot.hist2d(gas.x.value_in(units.AU),
                              gas.y.value_in(units.AU),
                              (200, 200), cmap=plt.cm.jet)
                pyplot.scatter(gas.x.value_in(units.AU), gas.y.value_in(units.AU), c='y', s=100)
                if False:
                  for gi in gas:
                    pyplot.arrow(gi.x.value_in(units.AU), gi.y.value_in(units.AU),
                             gi.vx.value_in(units.AU/units.yr), gi.vy.value_in(units.AU/units.yr),
                             head_width=0.05, head_length=0.1)

#                sph_particles_plot(gas, u_range=[min(gas.u), max(gas.u)], width=lim, alpha = 1)
#                if len(sinks):
#                    scatter(sinks.x.value_in(units.AU), sinks.y.value_in(units.AU), c='y', s=100)
                if image_id == snapshot_id:
                    break
    pyplot.xlabel("X [pc]")
    pyplot.ylabel("Y [pc]")
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "DiskWind.h5",
                      help="output filename [DiskWind.h5]")
    result.add_option("-l", unit=units.parsec,
                      dest="lim", type="float", default = None,
                      help="boxsize [%default]")
    result.add_option("-i", 
                      dest="image_id", type="int", default = -1,
                      help="image id [%default]")
    result.add_option("-p", 
                      dest="nice_plot", type="int", default = 1,
                      help="nice plot [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    filename = "hydro_GMC.h5"
    gas = read_set_from_file(filename, "hdf5")

    main(**o.__dict__)


