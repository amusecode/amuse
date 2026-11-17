import logging
import argparse

import numpy as np

from amuse.units import units, nbody_system
from amuse.datamodel import Particles

from amuse_metisse import Metisse

from amuse.ic.kroupa import new_kroupa_mass_distribution
# from amuse.ic.plummer import new_plummer_model

from amuse.io import write_set_to_file

logger = logging.getLogger("amuse")
# logger.setLevel(logging.DEBUG)
# logging.basicConfig(level=logging.DEBUG)


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument("--seed", type=int, default=127)
    result.add_argument("-n", "--number_of_stars", type=int, default=10000)
    result.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=100.0 | units.Myr,
        help="The time to evolve to",
    )
    result.add_argument(
        "-s",
        "--number_of_steps",
        type=int,
        default=100,
        help="Number of time steps (in logspace)",
    )
    result.add_argument(
        "-m",
        "--mass_min",
        type=units.MSun,
        default=0.1 | units.MSun,
        help="The minimum mass",
    )
    result.add_argument(
        "-M",
        "--mass_max",
        type=units.MSun,
        default=100.0 | units.MSun,
        help="The maximum mass",
    )
    result.add_argument(
        "--metallicity_dir",
        type=str,
        default="/Users/rieder/Code/UvA/Toonen/tres3.0/amuse/src/amuse_metisse/data/Hydrogen",
        help="The metallicity directory",
    )
    result.add_argument(
        "--metallicity_dir_he",
        type=str,
        default="/Users/rieder/Code/UvA/Toonen/tres3.0/amuse/src/amuse_metisse/data/Helium",
        help="The metallicity-He directory",
    )
    return result


def setup_metisse(metallicity_dir=None, metallicity_dir_he=None, **kwargs):
    # instance = Metisse(redirection="file")
    # instance = Metisse(redirection="none")
    instance = Metisse()

    # print(f"Setting parameters: {metallicity_dir=}, {metallicity_dir_he=}")
    # instance.parameters.metallicity_dir = metallicity_dir
    # instance.parameters.metallicity_dir_he = metallicity_dir_he

    # instance.parameters.wd_mass_scheme = "Modified_mestel"
    # instance.parameters.bhns_mass_scheme = "Belczynski2008"
    # instance.parameters.metallicity = 0.02

    # instance.commit_parameters()
    return instance


def test_metisse_sun(**kwargs):
    instance = setup_metisse(**kwargs)
    star = Particles(1)
    # star.mass = 4.7893208794726441 | units.MSun
    star.mass = 1.0 | units.MSun

    stars_in_metisse = instance.particles.add_particles(star)
    print(instance.parameters)
    assert stars_in_metisse[0].mass == 1.0 | units.MSun
    # assert stars_in_metisse[0].mass == 4.7893208794726441 | units.MSun
    print(stars_in_metisse[0])
    print("Evolving...")
    instance.evolve_model(1000.0 | units.yr)
    # instance.evolve_one_step(1)
    print(stars_in_metisse[0])
    print("Done")
    instance.stop()


def test_metisse_twostars(**kwargs):
    instance = setup_metisse(**kwargs)
    star = Particles(2)
    star.mass = [0.3, 2.5] | units.MSun

    stars_in_metisse = instance.particles.add_particles(star)
    print(stars_in_metisse[0])
    print("Evolving...")
    instance.evolve_model(1000.0 | units.yr)
    print(stars_in_metisse.stellar_type)
    print("Done")
    instance.stop()

def test_metisse_kroupa(**kwargs):
    instance = setup_metisse(**kwargs)
    number_of_stars = 1000
    star = Particles(number_of_stars)
    star.mass = new_kroupa_mass_distribution(number_of_stars, mass_min=0.3 | units.MSun, mass_max=10.0 | units.MSun)

    stars_in_metisse = instance.particles.add_particles(star)
    print(stars_in_metisse[0])
    print("Evolving...")
    instance.evolve_model(1000.0 | units.yr)
    print(stars_in_metisse)
    print("Done")
    instance.stop()


def evolve_stars_with_metisse(stars, age, **kwargs):
    instance = setup_metisse(**kwargs)
    stars_in_metisse = instance.particles.add_particles(stars)
    instance.evolve_model(age)
    stars_after_evolution = stars_in_metisse.copy()
    instance.stop()
    return stars_after_evolution


def evolve_stars_with_metisse_channel(stars, age, **kwargs):
    instance = setup_metisse(**kwargs)
    stars.mass = stars.mass_initial  # !!! PAY ATTENTION HERE !!!
    stars_in_metisse = instance.particles.add_particles(stars)
    if age > 0 | units.Myr:
        instance.evolve_model(age)
    print(stars[0])
    stars_in_metisse.new_channel_to(stars).copy_attributes(
        [
            "mass",
            "luminosity",
            "temperature",
            "radius",
        ]
    )
    print(stars[0])
    # import matplotlib.pyplot as plt
    # plt.scatter(
    #     stars.mass.value_in(units.MSun),
    #     stars.stellar_type.value_in(units.stellar_type),
    # )
    # ax = plt.gca()
    # ax.set_xscale("log")
    # #ax.set_yscale("log")
    # plt.show()
    # stars_after_evolution = stars_in_metisse.copy()
    instance.stop()
    return stars


def evolve_stars_metisse(
    number_of_stars=100,
    time_end=100.0 | units.Myr,
    number_of_steps=100,
    mass_min=0.1 | units.MSun,
    mass_max=100.0 | units.MSun,
    **kwargs
):
    mass = np.logspace(
        np.log10(mass_min.value_in(units.MSun)),
        np.log10(mass_max.value_in(units.MSun)),
        number_of_stars
    ) | units.MSun
    stars = Particles(number_of_stars)
    stars.mass = mass

    times = np.logspace(
        np.log10((time_end / (number_of_steps**1.5)).value_in(units.Myr)),
        np.log10(time_end.value_in(units.Myr)),
        number_of_steps
    ) | units.Myr
    i = 0
    for  time in times:
        print(f"Evolving to time: {time}")
        stars_metisse = evolve_stars_with_metisse(stars, time, **kwargs)
        stars_metisse.age = time
        write_set_to_file(
            stars_metisse, f"stars_metisse_{i:06d}.amuse"
        )
        i += 1

def evolve_stars_metisse_channels(
    number_of_stars=100,
    time_end=100.0 | units.Myr,
    number_of_steps=100,
    mass_min=0.1 | units.MSun,
    mass_max=100.0 | units.MSun,
    **kwargs
):
    mass = np.logspace(
        np.log10(mass_min.value_in(units.MSun)),
        np.log10(mass_max.value_in(units.MSun)),
        number_of_stars
    ) | units.MSun
    stars = Particles(number_of_stars)
    stars.mass_initial = mass

    times = np.logspace(
        np.log10((time_end / (number_of_steps**1.5)).value_in(units.Myr)),
        np.log10(time_end.value_in(units.Myr)),
        number_of_steps
    ) | units.Myr
    i = 0
    stars = evolve_stars_with_metisse_channel(stars, 0 | units.Myr, **kwargs)
    for  time in times:
        print(f"Evolving to time: {time}")
        stars = evolve_stars_with_metisse_channel(stars, time, **kwargs)
        # stars_metisse.age = time
        write_set_to_file(
            stars, f"stars3_metisse_{i:06d}.amuse"
        )
        i += 1

def evolve_stars_metisse_continuous(
    number_of_stars=100,
    time_end=100.0 | units.Myr,
    number_of_steps=100,
    mass_min=0.1 | units.MSun,
    mass_max=100.0 | units.MSun,
    **kwargs
):
    """
    Set up METISSE once and keep evolving the stars
    """
    mass = np.logspace(
        np.log10(mass_min.value_in(units.MSun)),
        np.log10(mass_max.value_in(units.MSun)),
        number_of_stars
    ) | units.MSun
    stars = Particles(number_of_stars)
    stars.mass_initial = mass
    stars.mass = mass

    times = np.logspace(
        np.log10((time_end / (number_of_steps**1.5)).value_in(units.Myr)),
        np.log10(time_end.value_in(units.Myr)),
        number_of_steps
    ) | units.Myr
    i = 0
    evo = setup_metisse(**kwargs)
    stars_in_metisse = evo.particles.add_particles(stars)
    channel_from_evo = stars_in_metisse.new_channel_to(stars)
    channel_from_evo.copy()
    for time in times:
        print(f"Evolving to time: {time}")
        evo.evolve_model(time)
        channel_from_evo.copy()
        # stars_metisse.age = time
        write_set_to_file(
            stars, f"stars6_metisse_{i:06d}.amuse"
        )
        i += 1


def main():
    args = new_argument_parser().parse_args()
    np.random.seed(args.seed)

    # test_metisse_sun(**vars(args))
    # test_metisse_twostars(**vars(args))
    # test_metisse_kroupa(**vars(args))
    # evolve_stars_metisse(start=0, **vars(args))
    # evolve_stars_metisse_channels(start=0, **vars(args))
    evolve_stars_metisse_continuous(start=0, **vars(args))

    # TODO:
    # METISSE tracks have a minimum mass, below that mass channels return zeros
    # So need to be able to return what this minimum (and maximum) mass is.


if __name__ == "__main__":
    main()
