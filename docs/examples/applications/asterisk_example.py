import numpy.random

from amuse.units import units, nbody_system
from amuse.ic.flatimf import new_flat_mass_distribution
from amuse.ic.plummer import new_plummer_model
from amuse.ext.particles_with_color import new_particles_with_blackbody_color
from amuse.community.seba.interface import SeBa
from amuse.community.bhtree.interface import BHTree
from amuse.community.asterisk.interface import Asterisk


def new_stellar_evolution(particles):
    stellar_evolution = SeBa()
    stellar_evolution.particles.add_particles(particles)
    return stellar_evolution


def new_gravity(particles, converter):
    gravity = BHTree(converter)
    gravity.particles.add_particles(particles)
    return gravity


if __name__ in ('__main__', '__plot__'):
    number_of_particles = 100

    # create a Plummer sphere with a number of stars
    numpy.random.seed(12345)
    masses = new_flat_mass_distribution(number_of_particles)
    converter = nbody_system.nbody_to_si(1.0 | units.parsec, masses.sum())
    particles = new_plummer_model(number_of_particles, converter)
    particles.mass = masses
    particles.move_to_center()

    # create simulation codes
    gravity = new_gravity(particles, converter)
    stellar_evolution = new_stellar_evolution(particles)

    # create channels to and from the local particle set and the simulations
    from_gravity_to_local = gravity.particles.new_channel_to(particles)
    from_stellar_evolution_to_local = \
        stellar_evolution.particles.new_channel_to(particles)
    from_stellar_evolution_to_local.copy()

    # creating colored particles
    particles = new_particles_with_blackbody_color(particles)
    particles.alpha = 1.0
    particles.radius = stellar_evolution.particles.radius.sqrt() * \
        (1e4 | units.parsec).sqrt()

    # creating visualization code
    converter = nbody_system.nbody_to_si(10.0 | units.parsec, masses.sum())
    visualization = Asterisk(converter, redirection="none")

    # optional: change OpenGL perspective settings
    # visualization.set_field_of_view(45.0)
    # visualization.set_z_near(0.1 | units.parsec)
    # visualization.set_z_far(3000.0 | units.parsec)

    # initialize code (creates gui windows)
    visualization.initialize_code()

    # optional: set the zoom and rotation of the visualization
    # visualization.parameters.rotation = (15, -15, 45)
    # visualization.parameters.camera_distance = 100 | units.parsec

    # add (now colored) particles to visualization
    visualization.particles.add_particles(particles)
    from_local_to_viz = particles.new_channel_to(visualization.particles)
    visualization.store_view(0 | units.Myr)

    # evolve module for some time
    for i in range(1, 100):
        print 'starting evolve to time = ', (i * 0.1 | units.Myr)
        target_time = i * 0.1 | units.Myr
        gravity.evolve_model(target_time)
        from_gravity_to_local.copy()
        stellar_evolution.evolve_model(target_time)
        from_stellar_evolution_to_local.copy()
        from_local_to_viz.copy_attributes(
            ["x", "y", "z", "red", "green", "blue"])
        visualization.particles.radius = (
                stellar_evolution.particles.radius.sqrt()
                * (1e4 | units.parsec).sqrt()
                )

        print 'updating visualization to time = ', target_time
        visualization.store_view(target_time)

    visualization.stop()
    gravity.stop()
    stellar_evolution.stop()
