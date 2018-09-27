"""
    Evolve the orbital evolution in the galactic center potential
    with different N-body codes.
    All codes use default settings and are integrated using Bridge.
    Since it is only a single particle being evolved, and the potential
    is handled by bridge (externally) there should not be any difference
    between the codes.
"""
from __future__ import print_function

import numpy
from itertools import cycle

from amuse.units import units, quantities, nbody_system
from amuse.units.optparse import OptionParser
from amuse.datamodel import Particle, Particles

from amuse.support.exceptions import AmuseException
from amuse.support.console import set_printing_strategy

from matplotlib import pyplot
from amuse import plot as aplot

from amuse.ext import static_potentials

from amuse.couple.bridge import Bridge
from amuse.community.ph4.interface import ph4
from amuse.community.huayno.interface import Huayno
from amuse.community.adaptb.interface import Adaptb
from amuse.community.bhtree.interface import BHTree
from amuse.community.tupan.interface import Tupan

codelist = {"huayno": Huayno, "ph4": ph4, "adaptb": Adaptb, "bhtree": BHTree,
            "tupan": Tupan}

coordinate_correction = [-8.08, 0, 6.68] | units.parsec


def subplot(potential, orbits, codes, fit_orbit, labels):
    hor, vert = labels
    X = numpy.linspace(-160, 160, 100) | units.parsec
    Y = numpy.linspace(-160, 160, 100) | units.parsec
    X, Y = quantities.meshgrid(X, Y)

    if labels == 'xy':
        pot_args = [X, Y, 0 | units.parsec]
        fit_horizontal = fit_orbit[0]
        fit_vertical = fit_orbit[1]
    elif labels == 'xz':
        pot_args = [X, 0 | units.parsec, Y]
        fit_horizontal = fit_orbit[0]
        fit_vertical = fit_orbit[2]
    elif labels == 'yz':
        pot_args = [0 | units.parsec, X, Y]
        fit_horizontal = fit_orbit[1]
        fit_vertical = fit_orbit[2]

    phi = potential.get_potential_at_point(None, *pot_args)
    aplot.imshow_color_plot(X, Y, phi, cmap="Blues")
    del aplot.UnitlessArgs.arg_units[2]

    aplot.scatter(0 | units.parsec, 0 | units.parsec, c='black')
    aplot.plot(fit_horizontal, fit_vertical, c="green", ls="--",
               label="Kruijssen et al. 2014")

    colors = cycle(['red', 'black', 'yellow', 'grey', 'green'])

    for full_orbit, code in zip(orbits, codes):
        horizontal = full_orbit.x if hor == 'x' else full_orbit.y
        vertical = full_orbit.y if vert == 'y' else full_orbit.z
        color = next(colors)
        aplot.plot(horizontal, vertical, c=color, label=code)
        aplot.scatter(horizontal[-1], vertical[-1], c=color, edgecolor=color)

    pyplot.axis('equal')
    aplot.xlim([-150, 150] | units.parsec)
    aplot.ylim([-150, 150] | units.parsec)
    aplot.xlabel(hor)
    aplot.ylabel(vert)


def plot_orbit(potential, orbits, codes, fit_orbit, filename):
    pyplot.figure(figsize=(10, 10))

    pyplot.subplot(221)
    subplot(potential, orbits, codes, fit_orbit, "xy")

    pyplot.subplot(222)
    subplot(potential, orbits, codes, fit_orbit, "xz")

    pyplot.subplot(223)
    subplot(potential, orbits, codes, fit_orbit, "yz")

    pyplot.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))

    if filename is None:
        pyplot.show()
    else:
        pyplot.savefig(filename)


def load_fit_orbit():
    (t, x, y, z, R, vx, vy, vz, vorb, l, b, vlos, mu_l, mub, mu_x,
     mu_y) = numpy.loadtxt("orbit_KDL15.dat", unpack=True)
    fit_orbit = [] | units.parsec
    for var in (x, y, z):
        fit_orbit.append(var | units.parsec)

    return fit_orbit


def setup_system():
    """
        Galacitic center potential and Arches cluster
        position from Kruijssen 2014
    """
    potential = static_potentials.Galactic_Center_Potential_Kruijssen()
    cluster = Particle()

    # At time 2.05 in KDL15
    cluster.position = [-17.55767, -53.26560, -9.39921] | units.parsec
    cluster.velocity = [-187.68008, 80.45276, 33.96556] | units.kms

    cluster.position += coordinate_correction
    return potential, cluster


def setup_bridge(potential, cluster, current_age, timestep, orbit, code):
    converter = nbody_system.nbody_to_si(
        current_age, cluster.position.length())
    gravity = codelist[code](converter)

    cluster.mass = 0 | units.MSun
    gravity.particles.add_particle(cluster)
    bridge = Bridge(timestep=timestep, use_threading=False)
    bridge.add_system(gravity, (potential,))

    gravity_to_orbit = gravity.particles.new_channel_to(orbit.particles)

    return bridge, gravity, gravity_to_orbit


def evolve_bridged_orbit_in_potential(potential, cluster, code, timestep,
                                      current_age):
    orbit = static_potentials.Position_In_Potential(potential, cluster)

    bridge, gravity, gravity_to_orbit = setup_bridge(
        potential, cluster, current_age, timestep, orbit, code)

    print(current_age, "->", orbit.particle.position)
    # Evolving backward
    gravity.particles.velocity *= -1
    while bridge.model_time < current_age-(timestep/2.):
        # print bridge.model_time, "-", gravity.particles.position[0]
        bridge.evolve_model(bridge.model_time + timestep)
        gravity_to_orbit.copy()
    gravity.particles.velocity *= -1

    gravity.model_time = 0 | units.Myr
    bridge.time = 0 | units.Myr
    cluster = gravity.particles[0]
    orbit.particle = cluster

    bridge, gravity, gravity_to_orbit = setup_bridge(
        potential, cluster, current_age, timestep, orbit, code)

    print(bridge.model_time, "->", orbit.particle.position)

    # Evolving forward
    full_orbit = Particles()
    while bridge.model_time < current_age-(timestep/2.):
        full_orbit.add_particle(orbit.particle.copy())
        bridge.evolve_model(bridge.model_time + timestep)
        gravity_to_orbit.copy()
    full_orbit.add_particle(orbit.particle.copy())
    print(bridge.model_time, "->", orbit.particle.position)

    full_orbit.position -= coordinate_correction
    return full_orbit


def evolve_orbit_and_plot(codes, filename, parameters):
    potential, cluster = setup_system()

    orbits = []
    for code in codes:
        print("evolving with", code)
        if code in codelist:
            full_orbit = evolve_bridged_orbit_in_potential(
                potential, cluster.copy(), code, **parameters)
        else:
            raise AmuseException("Unknow code: {}".format(code))
        orbits.append(full_orbit)
        print()
    fit_orbit = load_fit_orbit()

    plot_orbit(potential, orbits, codes, fit_orbit, filename)


def parse_arguments():
    parser = OptionParser()
    parser.add_option("-t", dest="timestep", type="float",
                      unit=units.Myr, default=0.1,
                      help="The timestep of the integrator [%default %unit].")
    parser.add_option("-a", dest="current_age", type="float",
                      unit=units.Myr, default=3.5,
                      help="The current age of the cluster [%default %unit].")
    parser.add_option("-c", dest="codes", action="append", type="string",
                      help="The code to evolve the particle (can be used"
                      "multiple times). Available codes: {}."
                      .format(codelist.keys()))
    parser.add_option("-f", dest="filename", type="string", default=None,
                      help="Save the plot in a file with 'name'"
                      "[show matplotlib if not given].")

    options, args = parser.parse_args()
    codes = options.codes or ["huayno"]
    filename = options.filename
    parameters = options.__dict__
    del parameters['codes']
    del parameters['filename']
    return codes, filename, parameters


if __name__ == "__main__":
    set_printing_strategy("custom", preferred_units=[units.MSun, units.parsec,
                                                     units.Myr, units.kms])
    codes, filename, parameters = parse_arguments()
    evolve_orbit_and_plot(codes, filename, parameters)
