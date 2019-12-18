import sys
import numpy
from numpy import random
from amuse.lab import *
from amuse.units.optparse import OptionParser
from amuse.units.quantities import as_vector_quantity
from matplotlib import pyplot

from prepare_figure import *
from distinct_colours import get_distinct

def movie_all(time, sun_and_planets):
    colors = ['r', 'b', 'g']

    print(time.in_(units.Gyr))
    R = [] | units.kpc
    for sp in sun_and_planets:
        R.append(sp.position.length())
    pyplot.subplot(2,2,1)
    pyplot.scatter(sun_and_planets.x.value_in(units.kpc), 
                   sun_and_planets.y.value_in(units.kpc), 
                   c=['k', 'r'], s=10, lw=0)
    pyplot.subplot(2,2,2)
    pyplot.scatter(R.value_in(units.kpc), 
                   sun_and_planets.z.value_in(units.kpc), 
                   c=['k', 'r'], s=10, lw=0)
    pyplot.xlabel("R [kpc]")
    pyplot.ylabel("Z [kpc]")

    R = [] | units.kpc
    R.append((sun_and_planets[1].position-sun_and_planets[0].position).length())
    pyplot.subplot(2,2,3)
    pyplot.scatter(-time.value_in(units.Gyr), 
                   R.value_in(units.kpc), 
                   c=['k', 'r'], s=10, lw=0)
    pyplot.xlabel("t [Gyr]")
    pyplot.ylabel("r [kpc]")

def new_option_parser():
    result = OptionParser()
    result.add_option("--seed", 
                      dest="seed", type="int", default = 666,
                      help="random number seed [%default]")
    result.add_option("-N", 
                      dest="Nstars", type="int", default = 3,
                      help="Number of stars [%default]")
    result.add_option("-f", 
                      dest="filename", default = "initial_cluster.amuse",
                      help="input filename")
    result.add_option("-R", unit=units.parsec,
                      dest="Rcluster", type="int", default = 1000|units.AU,
                      help="Cluster size [%default]")
    result.add_option("-n", 
                      dest="n", type="int", default = 100,
                      help="Number of pebbels per star [%default]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    filename = "sunandM67.hdf5"
    sun_and_planets = read_set_from_file(filename, "amuse")

    colors = get_distinct(4)
    figure = pyplot.figure(figsize=(16, 12))
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)

    x_label = "t [Gyr]"
    y_label = "$Z [kpc]$"
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    
#    pyplot.ylim(0, 1)
    R = [] 
    Z = [] | units.kpc
    d = [] | units.kpc
    t = [] | units.Gyr
    dmin = 100| units.kpc
    tmin = 0| units.Gyr
    zEmin = 0| units.kpc
    zMmin = 0| units.kpc
    rEmin = 0| units.kpc
    rMmin = 0| units.kpc
    for saps in sun_and_planets.history:
        R.append([saps[0].position.length().value_in(units.kpc), saps[1].position.length().value_in(units.kpc)])
        Z.append(saps.z)

        t.append(saps.get_timestamp())
        d.append((saps[1].position-saps[0].position).length())
    
        if d[-1]<dmin:
            dmin = d[-1]
            tmin = t[-1]
            zEmin = Z[-1][0]
            zMmin = Z[-1][0]
            rEmin = saps[0].position.length()
            rMmin = saps[1].position.length()
    print("minimal distance:", tmin.value_in(units.Myr), dmin.in_(units.parsec), rEmin.in_(units.parsec), rMmin.in_(units.parsec), zEmin.in_(units.parsec), zMmin.in_(units.parsec))

    plot_td = True
    plot_RZ = False
    if plot_RZ:
        pyplot.plot(R, Z.value_in(units.kpc), lw=2)
        pyplot.xlabel("R [kpc]")
        pyplot.ylabel("Z [kpc]")
        pyplot.xlim(8.0, 10.2)
        pyplot.ylim(-1.0, 1.0)
        pyplot.scatter(R[0], Z[0].value_in(units.kpc), lw=1, s=100, c=colors[0])
        pyplot.scatter(R[-1], Z[-1].value_in(units.kpc), lw=1, s=100, c=colors[1])
        pyplot.scatter(rMmin.value_in(units.kpc), zMmin.value_in(units.kpc), lw=1, s=100, c=colors[2])
    elif plot_td:
        pyplot.plot(-t.value_in(units.Gyr), d.value_in(units.kpc), lw=3, c=colors[0])
        pyplot.xlabel("t [Gyr]")
        pyplot.ylabel("d [kpc]")
        pyplot.ylim(0, 6)
        pyplot.xlim(-5, 0)
    pyplot.savefig("sun_and_M67")
    pyplot.show()

