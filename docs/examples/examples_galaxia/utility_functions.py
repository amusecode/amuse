# ---------------------------------------------------------------------------
# Script that shows the utility functions of galaxia.                       |
# These functions include: the computation of the components of the         |
# tidal tensor at a given point, the eigen-values and the tidal radius      |
# of a star cluster.                                                        |
# For a complete explanation of the possible parameters and models          |
# included in galaxia, we refer the reader to the file: user_manual_galaxia.|
# ---------------------------------------------------------------------------

from amuse.units import units
from amuse.community.galaxia.interface import BarAndSpirals3D
import numpy
from amuse.units.quantities import VectorQuantity


def galaxy_model(
        initial_phase_bar=(-20*numpy.pi)/180,
        initial_phase_spiral_arms=(-20*numpy.pi)/180,
        pattern_speed_spiral_arms=20 | (units.kms/units.kpc),
        amplitude_spiral_arms=1100 | (units.kms**2/units.kpc),
        number_of_arms=4,
        separation_locus_spiral_arms=3.12 | units.kpc,
        tangent_pitch_angle=0.227194425,
        pattern_speed_bar=50 | (units.kms/units.kpc),
        mass_bar=1.1e10 | units.MSun,
        semimajor_axis_bar=3.12 | units.kpc,
        axis_ratio_bar=0.37):
    # Model of the Galaxy.
    # In this example, the Galaxy has two-dimensional bar and spiral arms.
    # The spiral arms model is the TWA
    # The bar does not grow adiabatically
    # The axisymmetric component has its default values from Allen & Santillan
    # (1990).
    galaxy = BarAndSpirals3D()
    galaxy.parameters.bar_contribution = True
    galaxy.parameters.bar_phase = initial_phase_bar
    galaxy.parameters.omega_bar = pattern_speed_bar
    galaxy.parameters.mass_bar = mass_bar
    galaxy.parameters.aaxis_bar = semimajor_axis_bar
    galaxy.parameters.axis_ratio_bar = axis_ratio_bar
    galaxy.parameters.spiral_contribution = True
    galaxy.parameters.spiral_model = 0
    galaxy.parameters.spiral_phase = initial_phase_spiral_arms
    galaxy.parameters.omega_spiral = pattern_speed_spiral_arms
    galaxy.parameters.amplitude = amplitude_spiral_arms
    galaxy.parameters.rsp = separation_locus_spiral_arms
    galaxy.parameters.m = number_of_arms
    galaxy.parameters.tan_pitch_angle = tangent_pitch_angle
    galaxy.commit_parameters()
    return galaxy


def tidal_tensor(t, x, y, z, galaxy):
    Fxx, Fyx, Fzx, Fxy, Fyy, Fzy, Fxz, Fyz, Fzz = galaxy.get_tidal_tensor(
        t, x, y, z)
    return VectorQuantity.new_from_scalar_quantities(
            Fxx, Fyx, Fzx,
            Fxy, Fyy, Fzy,
            Fxz, Fyz, Fzz)


def eigen_values(t, x, y, z, galaxy):
    l1, l2, l3 = galaxy.get_eigen_values(t, x, y, z)
    return VectorQuantity.new_from_scalar_quantities(l1, l2, l3)


if __name__ in('__main__', '__plot__'):

    # Tidal tensor, eigen values and tidal radius of a star cluster
    # with 1000 MSun that is located at solar radius.
    # The local stellar density is also calculated.

    mass_cluster = 1000 | units.MSun
    time, x, y, z = (
            0 | units.Myr,
            8.5 | units.kpc,
            0 | units.kpc,
            0 | units.kpc,
            )
    MilkyWay = galaxy_model()

    Tij = tidal_tensor(time, x, y, z, MilkyWay)
    eigenValues = eigen_values(time, x, y, z, MilkyWay)
    tidal_radius = MilkyWay.get_tidal_radius(time, x, y, z, mass_cluster)
    local_stellar_density = MilkyWay.get_local_density(time, x, y, z)

    print(
            'tidal tensor:',
            Tij.value_in(100*units.kms**2/units.kpc**2),
            '\n', 'eigen values:',
            eigenValues.value_in(100*units.kms**2/units.kpc**2),
            '\n', 'tidal radius:',
            tidal_radius.value_in(units.parsec),
            '\n', 'local stellar density:',
            local_stellar_density.value_in(units.MSun/units.kpc**3)
            )
