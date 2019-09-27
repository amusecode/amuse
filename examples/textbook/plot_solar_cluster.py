"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
#from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from amuse.units.optparse import OptionParser
from time import sleep

from distinct_colours import get_distinct

def plot_single_image(particles, lim):
#    nullfmt   = NullFormatter()         # no labels
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.05
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    colors = get_distinct(4)

#    pyplot.rcParams.update({'font.size': 30})
    fig = pyplot.figure(figsize=(12,12))	
#    ax = pyplot.gca()
#    ax.minorticks_on() # switch on the minor ticks
#    ax.locator_params(nbins=3)

    time = particles.get_timestamp()
#    pyplot.title("Cluster at t="+str(time.in_(units.Gyr)))
    xy = pyplot.axes(rect_scatter)
    xy.text(11,11, "Galaxy at t="+str(time.in_(units.Gyr)),  ha='left', va='bottom')
    xz = pyplot.axes(rect_histx)
    yz = pyplot.axes(rect_histy)
    xy.set_xlabel("X [kpc]")
    xy.set_ylabel("Y [kpc]")
    xz.set_ylabel("Z [kpc]")
    yz.set_xlabel("Z [kpc]")

    xy.scatter([0.0], [0.0], s=200, c=colors[1], marker='+')
    xz.scatter([0.0], [0.0], s=200, c=colors[1], marker='+')
    yz.scatter([0.0], [0.0], s=200, c=colors[1], marker='+')
    xy.scatter([8.5], [0.0], s=100, c=colors[2], marker='o', lw=0)
    xz.scatter([8.5], [0.0], s=100, c=colors[2], marker='o', lw=0)
    yz.scatter([0.0], [0.0], s=100, c=colors[2], marker='o', lw=0)

#    axHistx.xaxis.set_major_formatter(nullfmt)
#    axHisty.yaxis.set_major_formatter(nullfmt)

    positions = particles.position
    x, y, z = positions.x.value_in(units.kpc), positions.y.value_in(units.kpc), positions.z.value_in(units.kpc)

    xy.scatter(x, y, c=colors[0], lw=0)
    xy.set_xlim( (-lim, lim) )
    xy.set_ylim( (-lim, lim) )
    xz.scatter(x, z, c=colors[0], lw=0)
    yz.scatter(z, y, c=colors[0], lw=0)
    xz.set_xlim( xy.get_xlim() )
#    yz.set_xlim( (-0.2*lim, 0.2*lim) )
#    yz.set_xlim( xy.get_xlim() )
    yz.set_ylim( xy.get_xlim() )
    yz.set_xlim( (-0.1*lim, 0.1*lim) )
    xz.set_ylim( (-0.1*lim, 0.1*lim) )
#    pyplot.show()
#    fig.savefig('test.png')
    fig.savefig('SolarSiblings_static_galaxy')
#    fig.savefig('test.eps')

def main(filename, lim=-1|units.parsec, image_id=-1):
    if image_id<0:
        pyplot.ion()
    mc = read_set_from_file(filename, "hdf5") # particles
    time = 0
    snapshot_id = 0
    for si in mc.history:
        snapshot_id += 1
        time = si.get_timestamp()
        print("Snapshot=", snapshot_id, time.in_(units.Gyr))
        if image_id<0 or image_id == snapshot_id:
            m = 1
            plot_single_image(si, lim.value_in(units.kpc))
        if image_id == snapshot_id:
            break
        if image_id<0:
            pyplot.draw()
            pyplot.cla()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "proto_solar_cluster.hdf5",
                      help="output filename [%default]")
    result.add_option("-l", unit=units.parsec, 
                      dest="lim", type="float", default = 10|units.kpc,
                      help="axis length [%default]")
    result.add_option("-i", 
                      dest="image_id", type="int", default = -1,
                      help="image id [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


