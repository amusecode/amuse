# -*- coding: utf-8 -*-
"""
Calculated theoretical pre-main-sequance evolutionary tracks for stars of
various masses.
After Iben, ApJ 141, 993, 1965
"""
from __future__ import print_function
from matplotlib import pyplot

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse import datamodel


def simulate_evolution_tracks(
        masses=[0.5, 1.0, 1.25, 1.5, 2.25, 3.0, 5.0, 9.0, 15.0] | units.MSun
):
    stellar_evolution = MESA()
    stellar_evolution.parameters.AGB_wind_scheme = 0
    stellar_evolution.parameters.RGB_wind_scheme = 0

    stars = datamodel.Particles(len(masses), mass=masses)
    stars = stellar_evolution.pre_ms_stars.add_particles(stars)

    data = {}

    for star in stars:
        luminosity = [] | units.LSun
        temperature = [] | units.K
        time = [] | units.yr

        print('Evolving pre main sequence star with')
        print(('    mass:', star.mass))
        print(('    luminosity:', star.luminosity))
        print(('    radius:', star.radius))

        while star.stellar_type == 17 | units.stellar_type:
            luminosity.append(star.luminosity)
            temperature.append(star.temperature)
            time.append(star.age)
            star.evolve_one_step()

        print(('Evolved pre main sequence star to:', star.stellar_type))
        print(('    age:', star.age))
        print(('    mass:', star.mass))
        print(('    luminosity:', star.luminosity))
        print(('    radius:', star.radius))
        print()

        stardata = {}
        stardata['luminosity'] = luminosity
        stardata['temperature'] = temperature
        stardata['time'] = time
        data[star.mass] = stardata

    return data


def plot_track(data):
    figure = pyplot.figure(figsize=(6, 8))
    plot = figure.add_subplot(1, 1, 1)
    plot.set_title('Hertzsprung-Russell diagram', fontsize=12)
    temp_unit = units.K
    luminosity_unit = units.LSun
    for mass, stardata in list(data.items()):
        temperature = stardata['temperature']
        luminosity = stardata['luminosity']
        plot.loglog(
            temperature[4:].value_in(temp_unit),
            luminosity[4:].value_in(luminosity_unit),
            marker="s")  # first few points show transient
        plot.text(
            1.25 * temperature[-1].value_in(temp_unit),
            0.5 * luminosity[-1].value_in(luminosity_unit),
            str(mass))
    plot.set_xlabel('Effective Temperature [' + str(temp_unit) + ']')
    plot.set_ylabel('Luminosity [' + str(luminosity_unit) + ']')
    plot.set_xlim(10**4.6, 10**3.5)
    plot.set_ylim(1.0e-2, 1.e5)
    pyplot.show()


if __name__ == "__main__":
    data = simulate_evolution_tracks(masses=[0.5, 1.0, 1.25] | units.MSun)
    plot_track(data)
