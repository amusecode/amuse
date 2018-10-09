# ---------------------------------------------------------------------------
# This script generates the rotation curve of the Milky Way using galaxia.  |
# For a complete explanation of the possible parameters and models          |
# included in galaxia, we refer the reader to the file: user_manual_galaxia.|
# ---------------------------------------------------------------------------

from amuse.units import units
from amuse.community.galaxia.interface import BarAndSpirals3D
import numpy
from matplotlib import pyplot


def galaxy_model_with_bar_and_spirals(
        mass_bar=1.4e10 | units.MSun,
        initial_phase_bar_and_SA=0.35,
        semimajor_axis_bar=3.12 | units.kpc,
        axis_ratio_bar=0.37,
        vertical_axis_bar=1 | units.kpc,
        amplitude_spirals=3e6 | (units.MSun/units.kpc**3),
        scale_length_spirals=2.6 | units.kpc,
        number_of_arms=2,
        tangent_pitch_angle=0.2773,
        scale_height_spirals=0.3 | units.kpc,
        fiducial_radius=8 | units.kpc):
    # Galaxy model that uses a three dimensional bar (no adiabatic growing)
    # and the C&G model for the spiral arms
    galaxy = BarAndSpirals3D()
    galaxy.parameters.bar_contribution = True
    galaxy.parameters.bar_phase = initial_phase_bar_and_SA
    galaxy.parameters.mass_bar = mass_bar
    galaxy.parameters.aaxis_bar = semimajor_axis_bar
    galaxy.parameters.axis_ratio_bar = axis_ratio_bar
    galaxy.parameters.caxis_bar = vertical_axis_bar
    galaxy.parameters.spiral_contribution = True
    galaxy.parameters.spiral_model = 1
    galaxy.parameters.spiral_phase = initial_phase_bar_and_SA
    galaxy.parameters.spiral_density_amplitude = amplitude_spirals
    galaxy.parameters.r_sigma = scale_length_spirals
    galaxy.parameters.m = number_of_arms
    galaxy.parameters.tan_pitch_angle = tangent_pitch_angle
    galaxy.parameters.scale_height = scale_height_spirals
    galaxy.parameters.fiducial_radius = fiducial_radius
    galaxy.commit_parameters()
    return galaxy


def galaxy_model_purely_axisymmetric():
    # Axisymmetric Galaxy model with default values
    galaxy = BarAndSpirals3D()
    galaxy.commit_parameters()
    return galaxy


def plot_rotation_curves(r, vc, vc1):
    figure = pyplot.figure(figsize=(6, 6))
    ax = figure.add_subplot(111)
    ax.plot(x.value_in(units.kpc), vc.value_in(units.kms),
            label='Model with bar and spiral arms')
    ax.plot(x.value_in(units.kpc), vc1.value_in(
        units.kms), label='Only axisymmetric')
    ax.set_xlabel('Galactocentric radius [Kpc]')
    ax.set_ylabel('Circular velocity [km/s]')
    ax.legend(loc='lower right')
    pyplot.show()


if __name__ in('__main__', '__plot__'):

    x, y, z = numpy.linspace(
        0.01, 15, 500) | units.kpc, 0 | units.kpc, 0 | units.kpc

    MilkyWay = galaxy_model_with_bar_and_spirals()
    circular_velocity = MilkyWay.get_velcirc(x, y, z)

    AxiGalaxy = galaxy_model_purely_axisymmetric()
    circular_velocity_axisymmetric_component = AxiGalaxy.get_velcirc(x, y, z)

    plot_rotation_curves(x, circular_velocity,
                         circular_velocity_axisymmetric_component)
