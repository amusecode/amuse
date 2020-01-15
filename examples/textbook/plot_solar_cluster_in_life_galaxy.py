import sys
import numpy
import matplotlib
from matplotlib import pyplot
from amuse.lab import *
from distinct_colours import get_distinct
from prepare_figure import single_frame

def  plot_galaxy_and_stars(galaxy, stars):
    
    colors = get_distinct(3)
    single_frame('X [kpc]', 'Y [kpc]')
    xlim = 10
    pyplot.xlim(-xlim, xlim)
    pyplot.ylim(-xlim, xlim)
    ax = pyplot.gca()

    import numpy as np
    import pandas as pd
    from scipy import stats, integrate
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(color_codes=True)

    lim = 10|units.kpc
    p = galaxy.select(lambda x: x<lim,["x"])
    p = p.select(lambda x: x>-lim,["x"])
    p = p.select(lambda y: y<lim,["y"])
    p = p.select(lambda y: y>-lim,["y"])
    p = p.select(lambda r: r.length()>5|units.kpc,["position"])
    x = p.x.value_in(units.kpc)
    y = p.y.value_in(units.kpc)
    sns.kdeplot(x, y, ax=ax, shade=True, n_levels=20, shade_lowest=False)
    m = 100*numpy.sqrt(stars.mass/stars.mass.max())
    pyplot.scatter(stars.x.value_in(units.kpc), stars.y.value_in(units.kpc), c=colors[0], s=m, lw=0)
    pyplot.savefig("SolarSiblings_life_galaxy")


def main(filename, tplot):
    sc = read_set_from_file(filename, "hdf5") # particles
    time = 0 | units.Myr
    snapshot_id = 0
    for si in sc.history:
        si = si.copy()
        snapshot_id += 1
        print(snapshot_id)
        if snapshot_id%2:
            print("convert")
            gc = si
        else:
            time = si.get_timestamp()
        if time>=tplot:
            break
        gc.move_to_center()
        print("Snapshot=", snapshot_id, time.in_(units.Gyr), len(gc), len(si))
    plot_galaxy_and_stars(gc, si)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "proto_solar_cluster_life.hdf5",
                      help="output filename [%default]")
    result.add_option("-t", unit=units.Gyr, 
                      dest="tplot", type="float", default = 3|units.Gyr,
                      help="moment of the plot [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


