import numpy as np
from amuse.datamodel import Particles
from amuse.units import nbody_system, units
from amuse.test.amusetest import TestWithMPI
from amuse.community.rebound.interface import ReboundInterface
from amuse.community.rebound.interface import Rebound


class ReboundInterfaceTests(TestWithMPI):

    def test_add_two_particles_serial_and_retrieve(self):
        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        instance.initialize_code()

        res1, error = instance.new_particle(
            mass=11.0, radius=2.0, x=1.0, y=2.0, z=3.0, vx=4.0, vy=5.0, vz=6.0
        )
        res2, error = instance.new_particle(
            mass=21.0, radius=5.0, x=10.0, y=0.0, z=0.0, vx=10.0, vy=0.0, vz=0.0
        )

        self.assertEqual(error, 0)
        self.assertEqual(res1, 0)
        self.assertEqual(res2, 1)

        m1, error = instance.get_mass(res1)
        self.assertEqual(error, 0)
        self.assertEqual(m1, 11)

        m2, error = instance.get_mass(res2)
        self.assertEqual(error, 0)
        self.assertEqual(m2, 21)

        x, y, z, error = instance.get_position(res1)
        self.assertEqual(error, 0)
        self.assertEqual(x, 1)
        self.assertEqual(y, 2)
        self.assertEqual(z, 3)

        vx, vy, vz, error = instance.get_velocity(res1)
        self.assertEqual(error, 0)
        self.assertEqual(vx, 4)
        self.assertEqual(vy, 5)
        self.assertEqual(vz, 6)
        instance.cleanup_code()
        instance.stop()

    def test_add_two_particles_parallel_and_retrieve(self):
        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        instance.initialize_code()

        instance.new_particle(
            [10, 20],
            [0, 0], [0, 0], [0, 0],
            [0, 0], [0, 0], [0, 0],
            [1, 1]
        )
        retrieved_state = instance.get_state(0)

        self.assertEqual(10.0,  retrieved_state['mass'])

        retrieved_state = instance.get_state([0, 1])
        self.assertEqual(10.0,  retrieved_state['mass'][0])
        self.assertEqual(20.0,  retrieved_state['mass'][1])
        instance.cleanup_code()
        instance.stop()

    def test_add_three_particles_and_retrieve(self):
        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        instance.initialize_code()

        indices, error = instance.new_particle(
            [10, 20, 30],
            [1, 2, 3], [0, 0, 0], [0, 0, 0],
            [0, 0, 0], [0, 0, 0], [0, 0, 0],
            [1, 0, 1]
        )
        n, error = instance.get_number_of_particles()
        self.assertEqual(n, 3)
        error = instance.delete_particle(indices[1])
        n, error = instance.get_number_of_particles()
        self.assertEqual(n, 2)
        instance.commit_particles()

        retrieved_state = instance.get_state([0, 2])
        self.assertEqual(10.0,  retrieved_state['mass'][0])
        self.assertEqual(30.0,  retrieved_state['mass'][1])
        mass, error = instance.get_mass(1)
        self.assertEqual(error, -1)
        instance.cleanup_code()
        instance.stop()

    def test_set_integrators(self):
        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        instance.initialize_code()
        error = 0
        list_of_integrators = [
            "ias15",
            "whfast",
            "sei",
            "leapfrog",
            # "hermes",  # removed
            # "whfast-helio",  # removed
            "none",
            "janus",
            "mercurius",
            "saba",
            "eos",
            "bs",
            "whfast512",
        ]
        for i, integrator in enumerate(list_of_integrators):
            print(i, integrator)
            retrieved_integrator = instance.get_integrator()
            self.assertEqual(error, 0)
            self.assertEqual(integrator, retrieved_integrator)
            try:
                error = instance.set_integrator(list_of_integrators[i+1])
                self.assertEqual(error, 0)
            except IndexError:
                retrieved_integrator = instance.get_integrator()
                self.assertEqual(error, 0)
                self.assertEqual(integrator, retrieved_integrator)

        instance.cleanup_code()
        instance.stop()

    def test_evolve_default_integrator(self):
        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        # Set up an equal-mass binary on a circular orbit:
        self.assertEqual(
            [0, 0],
            list(
                instance.new_particle(
                    0.5,  0.5, 0, 0,  0, 0.5, 0, 0.01
                ).values()
            )
        )
        self.assertEqual(
            [1, 0],
            list(
                instance.new_particle(
                    0.5,  -0.5, 0, 0,  0, -0.5, 0, 0.01
                ).values()
            )
        )
        self.assertEqual(0, instance.commit_particles())

        self.assertEqual(0, instance.evolve_model(np.pi))
        for result, expected in zip(
            list(instance.get_position(0).values()),
            [-0.5, 0.0, 0.0, 0]
        ):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(
            list(instance.get_position(1).values()),
            [0.5, 0.0, 0.0, 0]
        ):
            self.assertAlmostEqual(result, expected, 3)

        self.assertEqual(0, instance.evolve_model(2 * np.pi))
        for result, expected in zip(
            list(instance.get_position(0).values()),
            [0.5, 0.0, 0.0, 0]
        ):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(
            list(instance.get_position(1).values()),
            [-0.5, 0.0, 0.0, 0]
        ):
            self.assertAlmostEqual(result, expected, 3)

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test_fail_to_add_particle(self):
        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        index, err = instance.new_particle(
            0.5,
            0.5, 0, 0,
            0, 0.5, 0,
            0.01, 1
        )
        self.assertEqual(-10, err)

        instance.cleanup_code()
        instance.stop()

    def test_set_epsilon_squared(self):

        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        instance.initialize_code()

        instance.set_eps2(0.1 * 0.1)

        eps2 = instance.get_eps2()['epsilon_squared']

        self.assertEqual(0.1 * 0.1, eps2)
        instance.cleanup_code()
        instance.stop()

    def test_set_boundaries(self):

        instance = self.new_instance_of_an_optional_code(ReboundInterface)
        instance.initialize_code()

        instance.set_boundary("none")
        boundary_type = instance.get_boundary()
        self.assertEqual("none", boundary_type)

        instance.set_boundary("periodic")
        boundary_type = instance.get_boundary()
        self.assertEqual("periodic", boundary_type)

        instance.set_boundary_size(42.1)
        boundary_size = instance.get_boundary_size()['boundary_size']
        self.assertEqual(42.1, boundary_size)

        instance.cleanup_code()
        instance.stop()


class TestRebound(TestWithMPI):
    def evolve_earth_orbit(
        self,
        integrator=None,
    ):
        """
        Evolve Earth's orbit, with chosen integrator and other settings
        """
        convert_nbody = nbody_system.nbody_to_si(
            1.0 | units.MSun,
            149.5e6 | units.km,
        )

        interface = self.new_instance_of_an_optional_code(
            Rebound, convert_nbody
        )
        interface.initialize_code(redirection="none")
        interface.parameters.epsilon_squared = 0.0 | units.au**2
        if integrator is not None:
            interface.parameters.integrator = integrator
            self.assertEqual(
                interface.parameters.integrator,
                integrator,
            )
        else:
            self.assertEqual(
                interface.parameters.integrator,
                "ias15",
            )
        interface.dt_dia = 5000

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        interface.particles.add_particles(stars)

        interface.evolve_model(365.25 | units.day)
        interface.particles.copy_values_of_all_attributes_to(stars)

        position_at_start = earth.position.value_in(units.au)[0]
        position_after_full_rotation = earth.position.value_in(units.au)[0]
        self.assertAlmostRelativeEqual(
            position_at_start, position_after_full_rotation, 6
        )

        interface.evolve_model(365.25 + (365.25 / 2) | units.day)

        self.assertAlmostRelativeEqual(
            interface.model_time,
            365.25 + (365.25 / 2) | units.day, 3
        )
        interface.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.au)[0]
        self.assertAlmostRelativeEqual(
            -position_at_start, position_after_half_a_rotation, 3
        )

        interface.evolve_model(
            365.25 + (365.25 / 2) + (365.25 / 4) | units.day
        )
        self.assertAlmostRelativeEqual(
            interface.model_time,
            365.25 + (365.25 / 2) + (365.25 / 4) | units.day, 3
        )

        interface.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.au)[1]
        self.assertAlmostRelativeEqual(
            -position_at_start, position_after_half_a_rotation, 3
        )

        interface.cleanup_code()
        interface.stop()

    def new_system_of_sun_and_earth(self):
        stars = Particles(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(np.array((0.0, 0.0, 0.0)))
        sun.velocity = units.ms(np.array((0.0, 0.0, 0.0)))
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371)
        earth.position = units.km(np.array((149.5e6, 0.0, 0.0)))
        earth.velocity = units.ms(np.array((0.0, 29800, 0.0)))

        return stars

    def test_evolve_default_integrator(self):
        self.evolve_earth_orbit()

    def test_evolve_ias15_integrator(self):
        self.evolve_earth_orbit(integrator="ias15")

    def test_evolve_whfast_integrator(self):
        self.evolve_earth_orbit(integrator="whfast")

    def xtest_plot(self):
        # this test was meant to plot - but that's not a proper test
        convert_nbody = nbody_system.nbody_to_si(
            1.0 | units.MSun, 149.5e6 | units.km
        )

        instance = self.new_instance_of_an_optional_code(
            Rebound, convert_nbody
        )
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.au**2
        instance.dt_dia = 5000

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.particles.add_particles(stars)

        for x in range(1, 500, 10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_all_attributes_to(stars)
            stars.savepoint()

        instance.cleanup_code()
        instance.stop()

    def test_colliding_particles(self):
        particles = Particles(7)
        particles.mass = 0.001 | nbody_system.mass
        particles.radius = 0.1 | nbody_system.length
        particles.x = [
            -101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0
        ] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = (
            [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]]
        ) | nbody_system.speed

        instance = self.new_instance_of_an_optional_code(Rebound)
        instance.initialize_code()
        instance.parameters.set_defaults()
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        instance.evolve_model(1.0 | nbody_system.time)

        print(instance.particles)
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
        self.assertEqual(len(collisions.particles(0)), 3)
        self.assertEqual(len(collisions.particles(1)), 3)
        self.assertEqual(
            len(
                particles - collisions.particles(0) - collisions.particles(1)
            ),
            1
        )
        self.assertEqual(
            abs(collisions.particles(0).x - collisions.particles(1).x)
            <= (
                collisions.particles(0).radius
                + collisions.particles(1).radius
            ),
            [True, True, True]
        )

        sticky_merged = Particles(len(collisions.particles(0)))
        sticky_merged.mass = (
            collisions.particles(0).mass + collisions.particles(1).mass
        )
        sticky_merged.radius = collisions.particles(0).radius
        for p1, p2, merged in zip(
                collisions.particles(0),
                collisions.particles(1),
                sticky_merged
        ):
            merged.position = (p1 + p2).center_of_mass()
            merged.velocity = (p1 + p2).center_of_mass_velocity()

        print(instance.model_time)
        print(instance.particles)
        instance.particles.remove_particles(
            collisions.particles(0)
            + collisions.particles(1)
        )
        instance.particles.add_particles(sticky_merged)
        print(instance.particles)
        instance.evolve_model(1.0 | nbody_system.time)
        print()
        print(instance.model_time)
        print(instance.particles)
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 1.0 | nbody_system.time)
        self.assertEqual(len(collisions.particles(0)), 1)
        self.assertEqual(len(collisions.particles(1)), 1)
        self.assertEqual(
            len(
                instance.particles
                - collisions.particles(0)
                - collisions.particles(1)
            ), 2
        )
        self.assertEqual(
            abs(collisions.particles(0).x - collisions.particles(1).x)
            <= (
                collisions.particles(0).radius
                + collisions.particles(1).radius
            ),
            [True]
        )
        instance.stop()

    def test_energies(self):
        convert_nbody = nbody_system.nbody_to_si(
            1.0 | units.MSun, 149.5e6 | units.km
        )
        instance = self.new_instance_of_an_optional_code(
            Rebound, convert_nbody
        )
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.au**2
        instance.dt_dia = 5000

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.particles.add_particles(stars)
        self.assertAlmostRelativeEquals(
            instance.kinetic_energy, stars.kinetic_energy()
        )
        self.assertAlmostRelativeEquals(
            instance.potential_energy, stars.potential_energy(), 10
        )

        instance.stop()

    def test_out_of_box(self):
        instance = self.new_instance_of_an_optional_code(Rebound)
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = Particles(2)
        particles.position = ([0, 0, 0], [1, 0, 0]) | nbody_system.length
        particles.velocity = ([-1, 0, 0], [2, 0, 0]) | nbody_system.speed
        particles.radius = 0 | nbody_system.length
        particles.mass = 0.1 | nbody_system.mass

        instance.particles.add_particles(particles)
        instance.stopping_conditions.out_of_box_detection.enable()
        instance.parameters.stopping_conditions_out_of_box_size = 2 | nbody_system.length
        instance.parameters.stopping_conditions_out_of_box_use_center_of_mass = False
        instance.evolve_model(1 | nbody_system.time)
        self.assertTrue(
            instance.stopping_conditions.out_of_box_detection.is_set()
        )
        self.assertEqual(
            len(
                instance.stopping_conditions.out_of_box_detection.particles(0)
            ), 1
        )
        self.assertEqual(
            instance.stopping_conditions.out_of_box_detection.particles(0)[0].key,
            particles[1].key
        )
        instance.stop()

    def test_energies_nbody(self):
        instance = self.new_instance_of_an_optional_code(Rebound)
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = Particles(2)
        particles.position = ([0, 0, 0], [1, 0, 0]) | nbody_system.length
        particles.velocity = ([-1, 0, 0], [2, 0, 0]) | nbody_system.speed
        particles.radius = 0 | nbody_system.length
        particles.mass = 0.1 | nbody_system.mass

        instance.particles.add_particles(particles)
        self.assertAlmostRelativeEquals(
            instance.kinetic_energy, particles.kinetic_energy()
        )
        self.assertAlmostRelativeEquals(
            instance.potential_energy,
            particles.potential_energy(G=nbody_system.G)
        )
        instance.stop()

    def test_energies_subset(self):
        instance = self.new_instance_of_an_optional_code(Rebound)
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        subset0 = 0
        subset1 = instance.new_subset()
        print(subset1)
        self.assertEqual(subset1, 1)

        particles = Particles(2)
        particles.position = ([0, 0, 0], [1, 0, 0]) | nbody_system.length
        particles.velocity = ([-1, 0, 0], [2, 0, 0]) | nbody_system.speed
        particles.radius = 0 | nbody_system.length
        particles.mass = 0.1 | nbody_system.mass

        particles2 = Particles(2)
        particles2.position = ([0, 0, 0], [1, 0, 0]) | nbody_system.length
        particles2.velocity = ([-1, 0, 0], [2, 0, 0]) | nbody_system.speed
        particles2.radius = 0 | nbody_system.length
        particles2.mass = 2 | nbody_system.mass
        particles2.subset = subset1

        instance.particles.add_particles(particles)
        instance.particles.add_particles(particles2)
        print(instance.particles)
        self.assertAlmostRelativeEquals(
            instance.kinetic_energy, particles.kinetic_energy()
        )
        self.assertAlmostRelativeEquals(
            instance.potential_energy,
            particles.potential_energy(G=nbody_system.G)
        )
        self.assertAlmostRelativeEquals(
            instance.get_kinetic_energy(subset1),
            particles2.kinetic_energy()
        )
        self.assertAlmostRelativeEquals(
            instance.get_potential_energy(subset1),
            particles2.potential_energy(G=nbody_system.G)
        )
        # instance.stop()

    def test8(self):
        instance = self.new_instance_of_an_optional_code(Rebound)
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        subset0 = 0
        subset1 = instance.new_subset()
        print(subset1)
        self.assertEqual(subset1, 1)

        particles = Particles(2)
        particles.position = ([0, 0, 0], [1, 0, 0]) | nbody_system.length
        particles.velocity = ([0, 0, 0], [0, 0.5, 0]) | nbody_system.speed
        particles.radius = 0 | nbody_system.length
        particles.mass = 0.1 | nbody_system.mass

        particles2 = Particles(2)
        particles2.position = ([0, 0, 0], [2, 0, 0]) | nbody_system.length
        particles2.velocity = ([0, 0, 0], [0, 1, 0]) | nbody_system.speed
        particles2.radius = 0 | nbody_system.length
        particles2.mass = 1 | nbody_system.mass
        particles2.subset = subset1

        instance.particles.add_particles(particles)
        instance.particles.add_particles(particles2)
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostRelativeEquals(
            instance.kinetic_energy, particles.kinetic_energy(), 2
        )
        self.assertAlmostRelativeEquals(
            instance.potential_energy,
            particles.potential_energy(G=nbody_system.G),
            2
        )
        self.assertAlmostRelativeEquals(
            instance.get_kinetic_energy(subset1),
            particles2.kinetic_energy(),
            2
        )
        self.assertAlmostRelativeEquals(
            instance.get_potential_energy(subset1),
            particles2.potential_energy(G=nbody_system.G),
            2
        )
        particles_evolved = instance.particles.copy()
        instance.stop()

        instance1 = self.new_instance_of_an_optional_code(Rebound)
        instance1.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        particles1c = particles.copy()
        instance1.particles.add_particles(particles1c)
        instance1.evolve_model(0.1 | nbody_system.time)
        particles1evolved = particles_evolved[particles_evolved.subset == 0]
        print("HHH:", particles_evolved.subset)
        print("p2e:", particles1evolved)
        print(instance1.particles)
        self.assertAlmostRelativeEquals(
            instance1.particles.position,
            particles1evolved.position,
            10
        )
        self.assertAlmostRelativeEquals(
            instance1.particles.velocity,
            particles1evolved.velocity,
            10
        )
        instance1.stop()

        instance2 = self.new_instance_of_an_optional_code(Rebound)
        instance2.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        particles2c = particles2.copy()
        particles2c.subset = 0
        instance2.particles.add_particles(particles2c)
        instance2.evolve_model(0.1 | nbody_system.time)
        particles2evolved = particles_evolved[
            particles_evolved.subset == subset1
        ]
        print("HHH:", particles_evolved.subset)
        print("p2e:", particles2evolved)
        print(instance2.particles)
        self.assertAlmostRelativeEquals(
            instance2.particles.position,
            particles2evolved.position,
            10
        )
        self.assertAlmostRelativeEquals(
            instance2.particles.velocity,
            particles2evolved.velocity,
            10
        )
        instance2.stop()

    def test9(self):
        instance = self.new_instance_of_an_optional_code(Rebound)
        # instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = Particles(10)
        particles.x = np.arange(0, 10, 1) | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [0, 1, 0] | nbody_system.speed
        particles[5].velocity = [0, 0, 0] | nbody_system.speed
        particles.radius = 0 | nbody_system.length
        particles[5].mass = 1 | nbody_system.mass

        instance.particles.add_particles(particles)
        self.assertEqual(instance.particles[5].mass, 1 | nbody_system.mass)
        self.assertEqual(instance.particles[6].mass, 0 | nbody_system.mass)

        instance.evolve_model(1 | nbody_system.time)
        self.assertEqual(instance.particles[5].mass, 1 | nbody_system.mass)
        self.assertEqual(instance.particles[6].mass, 0 | nbody_system.mass)
        self.assertEqual(instance.particles[5].x, 5 | nbody_system.length)
        self.assertAlmostRelativeEquals(
            instance.particles[6].x, 5.5403023 | nbody_system.length, 5
        )
        instance.particles.remove_particle(particles[4])
        self.assertEqual(instance.particles[4].mass, 1 | nbody_system.mass)
        self.assertEqual(instance.particles[5].mass, 0 | nbody_system.mass)
        instance.evolve_model(1 | nbody_system.time)
        self.assertEqual(instance.particles[4].x, 5 | nbody_system.length)
        instance.stop()

    def test10(self):
        instance = self.new_instance_of_an_optional_code(Rebound)
        # instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = Particles(1)
        particles.x = 0 | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [0, 2, 0] | nbody_system.speed
        particles.radius = 1 | nbody_system.length
        particles.mass = 1 | nbody_system.mass

        instance.particles.add_particles(particles)
        instance.evolve_model(4 | nbody_system.time)
        self.assertAlmostRelativeEquals(
            instance.particles[0].velocity,
            [0, 2, 0] | nbody_system.speed
        )
        self.assertAlmostRelativeEquals(
            instance.particles[0].position,
            [0, 2 * 4, 0] | nbody_system.length,
            8
        )
        instance.stop()
