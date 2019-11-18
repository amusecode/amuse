# -*- coding: ascii -*-
from __future__ import print_function
from amuse.units import units, nbody_system
from amuse.datamodel import Particle
from amuse.community.athena.interface import Athena
from amuse.community.hermite.interface import Hermite

from matplotlib import pyplot


def hydro_grid_in_potential_well(mass=1 | units.MSun, length=100 | units.AU):
    converter = nbody_system.nbody_to_si(mass, length)

    # calculate density in field based on solar wind
    # gives a very low number
    molar_mass_hydrogen_proton = 1 | units.g / units.mol
    density_hydrogen_in_stellar_wind = 10 | 1 / units.cm**3
    particles_per_mol = 6.022e23 | 1 / units.mol
    density_hydrogen_in_stellar_wind_in_moles = (
        density_hydrogen_in_stellar_wind
        / particles_per_mol
    )
    density_gas = 100 * (
        density_hydrogen_in_stellar_wind_in_moles
        * molar_mass_hydrogen_proton
    ).as_quantity_in(units.MSun / units.AU**3)

    # override with higher number for plotting
    density_gas = 1e-3 | units.MSun / units.AU**3

    instance = Athena(converter)
    instance.initialize_code()
    instance.parameters.nx = 50
    instance.parameters.ny = 50
    instance.parameters.nz = 1
    instance.parameters.length_x = length
    instance.parameters.length_y = length
    instance.parameters.length_z = length
    instance.parameters.x_boundary_conditions = ("periodic", "periodic")
    instance.parameters.y_boundary_conditions = ("periodic", "periodic")
    instance.parameters.z_boundary_conditions = ("outflow", "outflow")

    # instance.stopping_conditions.number_of_steps_detection.enable()

    instance.set_has_external_gravitational_potential(1)

    instance.commit_parameters()

    grid_in_memory = instance.grid.copy()
    grid_in_memory.rho = density_gas
    pressure = 1 | units.Pa

    grid_in_memory.energy = pressure / (instance.parameters.gamma - 1)
    channel = grid_in_memory.new_channel_to(instance.grid)
    channel.copy()

    instance.initialize_grid()
    particle = Particle(
        mass=mass,
        position=length * [0.5, 0.5, 0.5],
        velocity=[0.0, 0.0, 0.0] | units.kms
    )

    gravity = Hermite(converter)
    dx = (grid_in_memory.x[1][0][0] - grid_in_memory.x[0]
          [0][0]).as_quantity_in(units.AU)
    gravity.parameters.epsilon_squared = dx**2
    gravity.particles.add_particle(particle)

    potential = gravity.get_potential_at_point(
        0 * instance.potential_grid.x.flatten(),
        instance.potential_grid.x.flatten(),
        instance.potential_grid.y.flatten(),
        instance.potential_grid.z.flatten()
    )

    potential = potential.reshape(instance.potential_grid.x.shape)
    instance.potential_grid.potential = potential

    instance.evolve_model(100 | units.yr)
    print(instance.get_timestep().value_in(units.yr))
    value_to_plot = instance.grid.rho[:, :, 0].value_in(
        units.MSun / units.AU**3)
    # value_to_plot = potential[...,...,0].value_in(potential.unit)
    plot_grid(value_to_plot)


def plot_grid(x):
    figure = pyplot.figure(figsize=(6, 6))
    plot = figure.add_subplot(1, 1, 1)
    mappable = plot.imshow(x, origin='lower')
    pyplot.colorbar(mappable)
    # figure.savefig('orszag_tang.png')
    pyplot.show()


if __name__ == '__main__':
    hydro_grid_in_potential_well()
