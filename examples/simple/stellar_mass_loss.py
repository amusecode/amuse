"""
    Set an artificial mass loss rate at a specific point
    in the stellar evolution using MESA
    and calculate the effect this has on the stellar radius.
"""
import numpy

from amuse.units import units
from amuse.datamodel import Particle, Particles
from amuse.support.console import set_printing_strategy

from amuse.community.mesa.interface import MESA


def print_report(star1, star2, mdot):
    report_string = ("mdot = {mdot}\n"
                     "  effective mdot = {eff_mdot}\n"
                     "  age: {star1.age} -> {star2.age}\n"
                     "    {d_age}\n"
                     "  radius: {star1.radius} -> {star2.radius}\n"
                     "    {d_radius}\n"
                     "  mass: {star1.mass} -> {star2.mass}\n"
                     "    {d_mass}\n"
                     "-> rdot: {rdot}\n\n")
    with open("radius_mass_loss_output.dat", 'a+') as f:
        f.write(report_string.format(
            star1=star1,
            star2=star2,
            mdot=mdot,
            eff_mdot=(star2.mass-star1.mass)/(star2.age-star1.age),
            rdot=(star2.radius-star1.radius)/(star2.age-star1.age),
            d_age=star2.age-star1.age,
            d_radius=star2.radius-star1.radius,
            d_mass=star2.mass-star1.mass))


def evolve_star_and_apply_mass_loss(radius, mdot):
    print("Evolve to radius = ", radius, "and then apply mdot =", mdot)

    stev = MESA(redirection='none')

    # We have to switch off all wind mass loss for manual mass loss to work
    stev.parameters.AGB_wind_scheme = 0
    stev.parameters.RGB_wind_scheme = 0

    star = stev.particles.add_particle(Particle(mass=2 | units.MSun))

    while star.radius < radius:
        star.evolve_one_step()
        print("evolved to:", star.age, "->", star.radius)

    star1 = star.copy()

    # High mass loss rates can only be calculated for small time steps
    star.time_step = 1. | units.yr
    star.mass_change = mdot
    print(star.mass_change)
    star.evolve_one_step()

    print_report(star1, star, mdot)


if __name__ == "__main__":
    set_printing_strategy(
            "custom",
            preferred_units=[
                units.RSun, units.MSun, units.Myr,
                units.MSun/units.yr,
                units.RSun/units.yr,
                ]
            )

    evolve_star_and_apply_mass_loss(
        39. | units.RSun, -1e-3 | units.MSun/units.yr)
