"""
    simple coupling of a stellar evolution and a gravitational code to simulate a cluster of stars.
"""

import os

from amuse.units.optparse import OptionParser
from amuse.units import units, nbody_system

from amuse.community.hermite0.interface import Hermite
from amuse.community.seba.interface import SeBa

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

from matplotlib import pyplot
from amuse import plot as aplot

from amuse.support.console import set_printing_strategy


def create_stars(number_of_stars, size):
    masses = new_salpeter_mass_distribution(number_of_stars, mass_min = 2|units.MSun)
    converter = nbody_system.nbody_to_si(masses.sum(), size)
    stars = new_plummer_model(number_of_stars, convert_nbody=converter)
    stars.mass = masses
    stars.zams_mass = masses

    return stars, converter

def plot_results(stars, time):
    mass_loss = stars.zams_mass - stars.mass

    x = stars.x.in_(units.parsec)
    y = stars.y.in_(units.parsec)

    pyplot.figure(figsize=(8,8))
    aplot.plot(x, y, "*")

    for x, y, mass_loss in zip(x.number, y.number, mass_loss):
        pyplot.annotate("%0.2f"%abs(mass_loss.number), xy=(x,y+2),
            horizontalalignment='center', verticalalignment='bottom')

    pyplot.axis('equal')
    pyplot.xlim([-60, 60])
    pyplot.ylim([-60, 60])
    aplot.xlabel("x")
    aplot.ylabel("y")
    pyplot.title("time = " + str(time))

    if not os.path.exists("plots"):
        os.mkdir("plots")

    name = "plots/plot_{0:=05}.png".format(int(time.value_in(units.Myr)))
    print "creating", name
    pyplot.savefig(name)
    pyplot.close()

def gravity_and_stellar_evolution(number_of_stars, size, end_time, sync_timestep=1|units.Myr, plot_timestep=10|units.Myr):

    stars, converter = create_stars(number_of_stars, size)

    gravity = Hermite(converter)
    stellar = SeBa()

    gravity.particles.add_particles(stars)
    stellar.particles.add_particles(stars)

    channels = []
    channels.append(stellar.particles.new_channel_to(gravity.particles, attributes=["mass", "radius"]))
    channels.append(stellar.particles.new_channel_to(stars))
    channels.append(gravity.particles.new_channel_to(stars))

    time = 0 | units.Myr
    while time <= end_time:
        for channel in channels:
            channel.copy()

        time += sync_timestep
        gravity.evolve_model(time)
        stellar.evolve_model(time)

        if time % plot_timestep < sync_timestep:
            plot_results(stars, time)

def parse_arguments():
    parser = OptionParser()
    parser.add_option("-N", dest="number_of_stars", type="int", default=100,
        help="The number of stars in the cluster [%default].")
    parser.add_option("-s", dest="size", type="float", unit=units.parsec, default=10,
        help="The total size of the cluster [%default %unit].")
    parser.add_option("-t", dest="end_time", type="float", unit=units.Gyr, default=0.1,
        help="The end time of the simulation [%default %unit].")

    options, args = parser.parse_args()
    return options.__dict__

if __name__ == "__main__":
    options = parse_arguments()
    set_printing_strategy("custom", preferred_units = [units.MSun, units.parsec, units.Myr], precision=3)
    gravity_and_stellar_evolution(**options)
