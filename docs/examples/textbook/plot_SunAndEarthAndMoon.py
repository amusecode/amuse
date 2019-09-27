import numpy
from matplotlib import pyplot
from amuse.lab import *

from prepare_figure import *
from distinct_colours import get_distinct

def orbital_elements(primary, secondary):
    stellar_pair = Particles()
    stellar_pair.add_particle(primary)
    stellar_pair.add_particle(secondary)
    Johannes.initialize_from_particles(stellar_pair)
    return Johannes.get_elements()

def read_and_process_file(filename):
    stars = read_set_from_file(filename, "hdf5")
    t = [] | units.yr
    Ek = [] | units.erg
    Ep = [] | units.erg
    for si in stars.history:
        Ek.append(si.kinetic_energy())
        Ep.append(si.potential_energy())
        time = si.get_timestamp()
        t.append(time)

    Etot = Ek + Ep
    ENorm = Etot[0]
    Etot = (Etot-ENorm)/ENorm
    return t, Etot
    
def main():
    colors = get_distinct(4)
    figure = pyplot.figure(figsize=(16, 12))
    ax = pyplot.gca()
    ax.minorticks_on() # switch on the minor ticks
    ax.locator_params(nbins=3)

    x_label = "t [year]"
    y_label = "$\Delta E/E_0$"
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.xlim(0, 110)

    t, E = read_and_process_file("SunAndEarthAndMoon_ISS.h5")
    pyplot.plot(t.value_in(units.yr), E, c=colors[0], label="ph4 ($dt_{param}=0.15$)")

    t, E = read_and_process_file("SunAndEarthAndMoon_TBwB.h5")
    pyplot.plot(t.value_in(units.yr), E, c=colors[3], label="Bridge with ph4 (E&M)")

    t, E = read_and_process_file("SunAndEarthAndMoon_TBB.h5")
    pyplot.plot(t.value_in(units.yr), E, c=colors[1], label="Three-body bridge")
    #    pyplot.plot(t.value_in(units.year), E)
    
    t, E = read_and_process_file("SunAndEarthAndMoon_TBBH.h5")
    pyplot.plot(t.value_in(units.yr), E, c=colors[2], label="Hierarchical bridge")

    pyplot.legend(loc=4, borderaxespad=0.)
    pyplot.savefig("fig_three_body_bridge")
    pyplot.show()

if __name__ == '__main__':
    main()












