import numpy as np
import pytest

from amuse.datamodel import Particles
from amuse.ext.orbital_elements import generate_binaries
from amuse.support.testing.amusetest import TestWithMPI
from amuse_tidymess.interface import Tidymess, TidymessInterface
from amuse.units import constants, nbody_system, units as u


class TestTidymessInterface(TestWithMPI):

    def test_initialization(self):
        """
        Test Tidymess initialization.
        """
        instance = self.new_instance_of_an_optional_code(TidymessInterface)
        assert instance is not None

        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test_setters_and_getters(self):
        """
        Test TidymessInterface setters and getters.
        """
        instance = self.new_instance_of_an_optional_code(TidymessInterface)
        assert instance is not None

        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())

        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 0)

        result = instance.new_particle(*np.arange(1, 16))

        self.assertEquals(result['index_of_the_particle'], 0.0)

        result = instance.get_state(0)
        self.assertEquals(result['mass'], 1.0)
        self.assertEquals(result['radius'], 2.0)
        self.assertEquals(result['x'], 3.0)
        self.assertEquals(result['y'], 4.0)
        self.assertEquals(result['z'], 5.0)
        self.assertEquals(result['vx'], 6.0)
        self.assertEquals(result['vy'], 7.0)
        self.assertEquals(result['vz'], 8.0)
        self.assertEquals(result['xi'], 9.0)
        self.assertEquals(result['kf'], 10.0)
        self.assertEquals(result['tau'], 11.0)
        self.assertEquals(result['wx'], 12.0)
        self.assertEquals(result['wy'], 13.0)
        self.assertEquals(result['wz'], 14.0)
        self.assertEquals(result['a_mb'], 15.0)

        instance.set_mass(0, 1.5)
        self.assertEquals(instance.get_mass(0)['mass'], 1.5)

        instance.set_position(0, 2.5, 3.5, 4.5)
        result = instance.get_position(0)
        self.assertEquals(result['x'], 2.5)
        self.assertEquals(result['y'], 3.5)
        self.assertEquals(result['z'], 4.5)

        instance.set_velocity(0, 5.5, 6.5, 7.5)
        result = instance.get_velocity(0)
        self.assertEquals(result['vx'], 5.5)
        self.assertEquals(result['vy'], 6.5)
        self.assertEquals(result['vz'], 7.5)

        instance.set_radius(0, 8.5)
        result = instance.get_radius(0)
        self.assertEquals(result['radius'], 8.5)

        instance.set_xi(0, 9.5)
        result = instance.get_xi(0)
        self.assertEquals(result['xi'], 9.5)

        instance.set_kf(0, 10.5)
        result = instance.get_kf(0)
        self.assertEquals(result['kf'], 10.5)

        instance.set_tau(0, 11.5)
        result = instance.get_tau(0)
        self.assertEquals(result['tau'], 11.5)

        instance.set_spin(0, 12.5, 13.5, 14.5)
        result = instance.get_spin(0)
        self.assertEquals(result['wx'], 12.5)
        self.assertEquals(result['wy'], 13.5)
        self.assertEquals(result['wz'], 14.5)

        instance.set_amb(0, 16)
        result = instance.get_amb(0)
        self.assertEquals(result['a_mb'], 16)

        instance.set_tidal_model(3)
        result = instance.get_tidal_model()
        self.assertEquals(result['tidal_model'], 3)

        instance.set_pn_order(2)
        result = instance.get_pn_order()
        self.assertEquals(result['pn_order'], 2)

        instance.set_magnetic_braking(1)
        result = instance.get_magnetic_braking()
        self.assertEquals(result['magnetic_braking'], 1)

        instance.set_speed_of_light(299792)
        result = instance.get_speed_of_light()
        self.assertEquals(result['speed_of_light'], 299792)

        instance.set_dt_mode(1)
        result = instance.get_dt_mode()
        self.assertEquals(result['dt_mode'], 1)

        instance.set_dt_const(0.025625)
        result = instance.get_dt_const()
        self.assertEquals(result['dt_const'], 0.025625)

        result = instance.get_time_step()
        self.assertEquals(result['time_step'], 0)

        instance.set_eta(0.625)
        result = instance.get_eta()
        self.assertEquals(result['eta'], 0.625)

        instance.set_n_iter(2)
        result = instance.get_n_iter()
        self.assertEquals(result['n_iter'], 2)

        instance.set_collision_mode(2)
        result = instance.get_collision_mode()
        self.assertEquals(result['collision_mode'], 2)

        instance.set_roche_mode(2)
        result = instance.get_roche_mode()
        self.assertEquals(result['roche_mode'], 2)

        instance.set_breakup_mode(1)
        result = instance.get_breakup_mode()
        self.assertEquals(result['breakup_mode'], 1)

        result = instance.get_num_integration_step()
        self.assertEquals(result['num_integration_step'], 0)

        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 1)

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test_adding_and_deleting_particles(self):
        """
        Test TidymessInterface creating and deleting particles.
        """
        instance = self.new_instance_of_an_optional_code(TidymessInterface)
        assert instance is not None

        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())

        # FIRST PARTICLE
        # --------------
        # initialize new particle with all attributes set to 1
        result = instance.new_particle(*np.ones(15)*1)
        self.assertEquals(result['index_of_the_particle'], 0)

        # SECOND PARTICLE
        # ---------------
        # initialize new particle with all attributes set to 1.1
        result = instance.new_particle(*np.ones(15)*1.1)
        self.assertEquals(result['index_of_the_particle'], 1)

        # check number of particles
        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 2)

        # check that indexes are correct
        first = instance.get_index_of_first_particle()
        self.assertEquals(first['index_of_the_particle'], 0)

        next = instance.get_index_of_next_particle(
            first['index_of_the_particle']
        )
        self.assertEquals(next['index_of_the_next_particle'], 1)

        # DELETE SECOND PARTICLE
        # ----------------------
        instance.delete_particle(1)

        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 1)

        first = instance.get_index_of_first_particle()
        self.assertEquals(first['index_of_the_particle'], 0)

        # THIRD PARTICLE
        # --------------
        result = instance.new_particle(*np.ones(15)*2)
        self.assertEquals(result['index_of_the_particle'], 2)

        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 2)

        first = instance.get_index_of_first_particle()
        self.assertEquals(first['index_of_the_particle'], 0)

        next = instance.get_index_of_next_particle(
            first['index_of_the_particle']
        )
        self.assertEquals(next['index_of_the_next_particle'], 2)

        result = instance.get_state(0)
        self.assertEquals(result['mass'], 1.0)
        self.assertEquals(result['x'], 1.0)
        self.assertEquals(result['y'], 1.0)
        self.assertEquals(result['z'], 1.0)
        self.assertEquals(result['vx'], 1.0)
        self.assertEquals(result['vy'], 1.0)
        self.assertEquals(result['vz'], 1.0)
        self.assertEquals(result['radius'], 1.0)
        self.assertEquals(result['xi'], 1.0)
        self.assertEquals(result['kf'], 1.0)
        self.assertEquals(result['tau'], 1.0)
        self.assertEquals(result['wx'], 1.0)
        self.assertEquals(result['wy'], 1.0)
        self.assertEquals(result['wz'], 1.0)
        self.assertEquals(result['a_mb'], 1.0)

        result = instance.get_state(2)
        self.assertEquals(result['mass'], 2.0)
        self.assertEquals(result['x'], 2.0)
        self.assertEquals(result['y'], 2.0)
        self.assertEquals(result['z'], 2.0)
        self.assertEquals(result['vx'], 2.0)
        self.assertEquals(result['vy'], 2.0)
        self.assertEquals(result['vz'], 2.0)
        self.assertEquals(result['radius'], 2.0)
        self.assertEquals(result['xi'], 2.0)
        self.assertEquals(result['kf'], 2.0)
        self.assertEquals(result['tau'], 2.0)
        self.assertEquals(result['wx'], 2.0)
        self.assertEquals(result['wy'], 2.0)
        self.assertEquals(result['wz'], 2.0)
        self.assertEquals(result['a_mb'], 2.0)

        instance.stop()

    def test_evolve_model(self):
        """
        Test TidymessInterface evolve_model with an equal mass binary.
        """
        instance = self.new_instance_of_an_optional_code(TidymessInterface)
        assert instance is not None

        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())

        self.assertEquals([0, 0], list(instance.new_particle(
            0.5, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0
        ).values()))
        self.assertEquals([1, 0], list(instance.new_particle(
            0.5, 0, -0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0
        ).values()))
        self.assertEquals(0, instance.commit_particles())

        self.assertEquals(0, instance.evolve_model(np.pi))  # half an orbit
        for result, expected in zip(instance.get_position(0).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 5)

        self.assertEquals(0, instance.evolve_model(2 * np.pi))  # full orbit
        for result, expected in zip(instance.get_position(0).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 5)

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()


class TestTidymess(TestWithMPI):

    def generate_HD80606b_system(self):
        """
        Initial conditions for the exoplanet system
        HD80606b. Initial conditions courtesy of
        Dr. Tjarda Boekholt.
        """
        HD80606 = Particles(2)
        star, planet = generate_binaries(
            primary_mass=2.0088092e30 | u.kg,
            secondary_mass=7.7459434e27 | u.kg,
            semi_major_axis=0.455 | u.au,
            eccentricity=0.9330
        )
        HD80606[0].name = 'HD80606'
        HD80606[0].mass = star.mass
        HD80606[0].radius = 702455.0 | u.km
        HD80606[0].x = star.x.as_quantity_in(u.km)
        HD80606[0].y = star.y.as_quantity_in(u.km)
        HD80606[0].z = star.z.as_quantity_in(u.km)
        HD80606[0].vx = star.vx.as_quantity_in(u.kms)
        HD80606[0].vy = star.vy.as_quantity_in(u.kms)
        HD80606[0].vz = star.vz.as_quantity_in(u.kms)
        HD80606[0].xi = 0.07
        HD80606[0].kf = 0.0
        HD80606[0].tau = 0.0 | u.yr
        HD80606[0].wx = 0.0 | 1 / u.s
        HD80606[0].wy = 0.0 | 1 / u.s
        HD80606[0].wz = 2.97188607137e-6 | 1 / u.s
        HD80606[0].a_mb = 0 | u.s

        HD80606[1].name = 'HD80606b'
        HD80606[1].mass = planet.mass
        HD80606[1].radius = 68488.3446 | u.km
        HD80606[1].x = planet.x.as_quantity_in(u.km)
        HD80606[1].y = planet.y.as_quantity_in(u.km)
        HD80606[1].z = planet.z.as_quantity_in(u.km)
        HD80606[1].vx = planet.vx.as_quantity_in(u.kms)
        HD80606[1].vy = planet.vy.as_quantity_in(u.kms)
        HD80606[1].vz = planet.vz.as_quantity_in(u.kms)
        HD80606[1].xi = 0.25
        HD80606[1].kf = 0.5
        HD80606[1].tau = 1e-2 | u.yr
        HD80606[1].wx = 0.0 | 1 / u.s
        HD80606[1].wy = 0.0 | 1 / u.s
        HD80606[1].wz = 0.000145444104333 | 1 / u.s
        HD80606[1].a_mb = 0 | u.s

        return HD80606

    def generate_figure8_system(self):
        """
        Generate initial conditions for a triple system
        in a figure 8 configuration.

        This is one of the solutions to the 3-body problem.
        The units are in dimensionless N-Body units.
        Initial conditions courtesy of Dr. Tjarda Boekholt.
        """
        figure8 = Particles(3)
        figure8[0].name = 'Star1'
        figure8[0].mass = 1 | nbody_system.mass
        figure8[0].radius = 5e-2 | nbody_system.length
        figure8[0].x = 0 | nbody_system.length
        figure8[0].y = 0 | nbody_system.length
        figure8[0].z = 0 | nbody_system.length
        figure8[0].vx = -0.93240737 | nbody_system.speed
        figure8[0].vy = -0.86473146 | nbody_system.speed
        figure8[0].vz = 0 | nbody_system.speed
        figure8[0].xi = 0.07
        figure8[0].kf = 0.02
        figure8[0].tau = 1e-2 | nbody_system.time
        figure8[0].wx = 0 | 1 / nbody_system.time
        figure8[0].wy = 0 | 1 / nbody_system.time
        figure8[0].wz = 0 | 1 / nbody_system.time
        figure8[0].a_mb = 0 | nbody_system.time

        figure8[1].name = 'Star2'
        figure8[1].mass = 1 | nbody_system.mass
        figure8[1].radius = 5e-2 | nbody_system.length
        figure8[1].x = 0.9700436 | nbody_system.length
        figure8[1].y = -0.24308753 | nbody_system.length
        figure8[1].z = 0 | nbody_system.length
        figure8[1].vx = 0.466203685 | nbody_system.speed
        figure8[1].vy = 0.43236573 | nbody_system.speed
        figure8[1].vz = 0 | nbody_system.speed
        figure8[1].xi = 0.07
        figure8[1].kf = 0.02
        figure8[1].tau = 1e-2 | nbody_system.time
        figure8[1].wx = 0 | 1 / nbody_system.time
        figure8[1].wy = 0 | 1 / nbody_system.time
        figure8[1].wz = 0 | 1 / nbody_system.time
        figure8[1].a_mb = 0 | nbody_system.time

        figure8[2].name = 'Star3'
        figure8[2].mass = 1 | nbody_system.mass
        figure8[2].radius = 5e-2 | nbody_system.length
        figure8[2].x = -0.9700436 | nbody_system.length
        figure8[2].y = 0.24308753 | nbody_system.length
        figure8[2].z = 0 | nbody_system.length
        figure8[2].vx = 0.466203685 | nbody_system.speed
        figure8[2].vy = 0.43236573 | nbody_system.speed
        figure8[2].vz = 0 | nbody_system.speed
        figure8[2].xi = 0.07
        figure8[2].kf = 0.02
        figure8[2].tau = 1e-2 | nbody_system.time
        figure8[2].wx = 0 | 1 / nbody_system.time
        figure8[2].wy = 0 | 1 / nbody_system.time
        figure8[2].wz = 0 | 1 / nbody_system.time
        figure8[2].a_mb = 0 | nbody_system.time

        return figure8

    def test_parameters_and_defaults(self):
        """
        Test Tidymess parameters attribute and their defaults.

        Although the default values are technically defined under
        parameters in the Tidymess interface.py, in practice these
        defaults are overridden by the defaults set in the Tidymess
        source code. The defaults in the interface.py were chosen to
        match the defaults set in the Tidymess standalone package.
        """
        system = self.generate_HD80606b_system()
        converter = nbody_system.nbody_to_si(
            system.mass.sum(), system[1].position.length()
        )
        instance = self.new_instance_of_an_optional_code(
            Tidymess, converter
        )
        assert instance is not None

        self.assertEquals(instance.parameters.tidal_model, 0)
        instance.parameters.tidal_model = 4
        self.assertEquals(instance.parameters.tidal_model, 4)

        self.assertEquals(instance.parameters.pn_order, 0)
        instance.parameters.pn_order = 1
        self.assertEquals(instance.parameters.pn_order, 1)

        self.assertEquals(instance.parameters.magnetic_braking, 0)
        instance.parameters.magnetic_braking = 1
        self.assertEquals(instance.parameters.magnetic_braking, 1)

        self.assertEquals(instance.parameters.speed_of_light, constants.c)
        instance.parameters.speed_of_light = 2e100 | u.m/u.s
        self.assertAlmostRelativeEquals(
            instance.parameters.speed_of_light.as_quantity_in(u.m / u.s),
            2e100 | u.m/u.s,
        )
        self.assertEquals(instance.parameters.dt_mode, 1)
        instance.parameters.dt_mode = 2
        self.assertEquals(instance.parameters.dt_mode, 2)

        dt_const = instance.unit_converter.to_si(0.015625 | nbody_system.time)
        self.assertAlmostEquals(instance.parameters.dt_const, dt_const, 3)
        instance.parameters.dt_const = 0.025625 | u.s
        self.assertAlmostEquals(instance.parameters.dt_const, 0.025625 | u.s, 3)

        self.assertEquals(instance.parameters.eta, 0.0625)
        instance.parameters.eta = 0.0726
        self.assertEquals(instance.parameters.eta, 0.0726)

        self.assertEquals(instance.parameters.n_iter, 1)
        instance.parameters.n_iter = 2
        self.assertEquals(instance.parameters.n_iter, 2)

        self.assertEquals(instance.parameters.collision_mode, 0)
        instance.parameters.collision_mode = 1
        self.assertEquals(instance.parameters.collision_mode, 1)

        self.assertEquals(instance.parameters.roche_mode, 0)
        instance.parameters.roche_mode = 2
        self.assertEquals(instance.parameters.roche_mode, 2)

        self.assertEquals(instance.parameters.breakup_mode, 0)
        instance.parameters.breakup_mode = 1
        self.assertEquals(instance.parameters.breakup_mode, 1)

        self.assertEquals(instance.parameters.initial_shape, 0)
        instance.parameters.initial_shape = 1
        self.assertEquals(instance.parameters.initial_shape, 1)

        instance.commit_parameters()

        instance.particles.add_particles(system)

        # check that parameters are still set
        # correctly after adding a particle
        self.assertEquals(instance.parameters.tidal_model, 4)
        self.assertEquals(instance.parameters.pn_order, 1)
        self.assertEquals(instance.parameters.magnetic_braking, 1)
        self.assertAlmostRelativeEquals(
            instance.parameters.speed_of_light.as_quantity_in(u.m / u.s),
            2e100 | u.m/u.s,
        )
        self.assertEquals(instance.parameters.dt_mode, 2)
        self.assertAlmostEquals(instance.parameters.dt_const, 0.025625 | u.s, 3)
        self.assertEquals(instance.parameters.eta, 0.0726)
        self.assertEquals(instance.parameters.n_iter, 2)
        self.assertEquals(instance.parameters.collision_mode, 1)
        self.assertEquals(instance.parameters.roche_mode, 2)
        self.assertEquals(instance.parameters.breakup_mode, 1)
        self.assertEquals(instance.parameters.initial_shape, 1)

        instance.stop()

    def test_adding_and_deleting_particles(self):
        """Test adding and deleting particles in Tidymess."""
        system1 = self.generate_HD80606b_system()
        converter = nbody_system.nbody_to_si(
            system1.mass.sum(), system1[1].position.length()
        )

        instance = self.new_instance_of_an_optional_code(Tidymess, converter)
        assert instance is not None

        instance.parameters.tidal_model = 0
        instance.commit_parameters()

        instance.particles.add_particles(system1)

        self.assertEquals(instance.get_number_of_particles(), 2)

        self.assertEquals(instance.model_time, 0 | u.s)
        self.assertAlmostRelativeEquals(
            instance.get_total_mass(),
            system1[0].mass + system1[1].mass,
            places=6
        )
        self.assertAlmostRelativeEquals(
            instance.get_total_radius(), system1[1].position.length()
        )

        instance.particles.remove_particle(system1[1])
        self.assertEquals(instance.get_number_of_particles(), 1)
        self.assertAlmostRelativeEquals(instance.get_total_mass(), system1[0].mass)

        system2 = self.generate_HD80606b_system()
        instance.particles.add_particles(system2)

        self.assertEquals(instance.get_number_of_particles(), 3)

        self.assertEquals(instance.particles[0], system1[0])
        self.assertEquals(instance.particles[1], system2[0])
        self.assertEquals(instance.particles[2], system2[1])

        instance.stop()

    def test_converting_spin_vectors(self):
        """
        Test the function for converting spin vectors.

        The comparison results were computed from an identical
        Tidymess standalone simulation.
        """
        system = self.generate_HD80606b_system()
        converter = nbody_system.nbody_to_si(
            system.mass.sum(), system[1].position.length()
        )
        instance = self.new_instance_of_an_optional_code(Tidymess, converter)
        assert instance is not None

        lod = 24.47 | u.day
        obl = 0 | u.deg
        psi = 0 | u.deg
        spin1 = instance.convert_spin_vectors_to_inertial(lod, obl, psi)
        self.assertEquals(spin1[0], 0 | 1/u.s)
        self.assertAlmostEquals(spin1[1], 0 | 1/u.s)
        self.assertAlmostEquals(spin1[2], 2.9718860713702659e-06 | 1/u.s)

        lod = 24.47 | u.day
        obl = 2.5 | u.deg
        psi = 3.0 | u.deg
        spin2 = instance.convert_spin_vectors_to_inertial(lod, obl, psi)
        self.assertAlmostRelativeEquals(
            spin2[0].value_in(1/u.yr), 2.1409959654769711e-1, places=3
        )
        self.assertAlmostRelativeEquals(
            spin2[1].value_in(1/u.yr), -4.0852636665140709, places=3
        )
        self.assertAlmostRelativeEquals(
            spin2[2].value_in(1/u.yr), 9.3696328783584647e+1, places=3
        )

        # check radians
        lod = 24.47 | u.day
        obl = (2.5 | u.deg).as_quantity_in(u.rad)
        psi = (3.0 | u.deg).as_quantity_in(u.rad)
        spin3 = instance.convert_spin_vectors_to_inertial(lod, obl, psi)
        self.assertAlmostRelativeEquals(
            spin3[0].value_in(1/u.yr), 2.1409959654769711e-1, places=3
        )
        self.assertAlmostRelativeEquals(
            spin3[1].value_in(1/u.yr), -4.0852636665140709, places=3
        )
        self.assertAlmostRelativeEquals(
            spin3[2].value_in(1/u.yr), 9.3696328783584647e+1, places=3
        )

        lod = 0.5 | u.day
        obl = 0 | u.deg
        psi = 0 | u.deg
        spin4 = instance.convert_spin_vectors_to_inertial(lod, obl, psi)
        self.assertEquals(spin4[0], 0 | 1/u.s)
        self.assertAlmostEquals(spin4[1], 0 | 1/u.s)
        self.assertAlmostEquals(spin4[2], 1.4544410433286076e-04 | 1/u.s)

        instance.stop()

    def test_begin_time(self):
        """
        Test that setting begin_time correctly creates a time offset.
        """
        dt = 5 | u.yr
        begin_time = (50 | u.yr).as_quantity_in(u.s)
        end_time = begin_time + dt

        system = self.generate_HD80606b_system()
        converter = nbody_system.nbody_to_si(
            system.mass.sum(), system[1].position.length()
        )
        instance = self.new_instance_of_an_optional_code(Tidymess, converter)
        assert instance is not None

        instance.parameters.tidal_model = 0
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.commit_parameters()

        self.assertEquals(instance.get_begin_time(), instance.model_time)
        self.assertEquals(instance.get_begin_time(), 0 | u.s)

        instance.set_begin_time(begin_time)
        instance.particles.add_particles(system)

        instance.evolve_model(instance.model_time + dt)

        self.assertAlmostEquals(instance.get_begin_time(), begin_time, places=3)
        self.assertAlmostEquals(instance.model_time, end_time, places=3)

        instance.stop()

    def test_evolving_backwards_in_time(self):
        """
        Test evolving backwards in time.

        2 identical figure 8 systems are created and evolved to t1.
        The second system is then evolved to t2 where t2 > t1.
        Finally the second system is evolved 'backwards' back to t1.
        System 2 should be identical to system 1 at the end of the
        backwards evolution.
        """
        t1 = 1e3 | nbody_system.time
        t2 = 2e3 | nbody_system.time
        self.assertFalse(t1 > t2)

        # System 1
        # --------
        system1 = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 1
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.commit_parameters()

        instance.particles.add_particles(system1)
        channel1 = instance.particles.new_channel_to(system1)

        instance.evolve_model(t1)
        channel1.copy()
        self.assertEquals(instance.model_time, t1)

        instance.stop()

        # System 2
        # --------
        system2 = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 1
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.commit_parameters()

        instance.particles.add_particles(system2)
        channel2 = instance.particles.new_channel_to(system2)

        instance.evolve_model(t2)
        channel2.copy()
        self.assertEquals(instance.model_time, t2)

        instance.evolve_model(t1)
        channel2.copy()
        self.assertEquals(instance.model_time, t1)

        instance.stop()

        attributes = [
            'x', 'y', 'z',
            'vx', 'vy', 'vz',
            'wx', 'wy', 'wz',
            'mass', 'radius',
            'xi', 'kf', 'tau'
        ]

        for attr in attributes:
            self.assertAlmostEquals(
                getattr(system1, attr),
                getattr(system2, attr)
            )

    def test_physical_units(self):
        """
        Evolve a system with tidal effects in physical units.

        The reference values used in this test were generated
        using the standalone Tidymess code. Because Tidymess
        applies slightly different unit conversion factors in
        physical-unit mode, the initial conditions were first
        converted to N-body units before running the simulation
        in Tidymess. The resulting outputs were then converted
        back to physical units for comparison with the AMUSE results.
        """
        system = self.generate_HD80606b_system()
        converter = nbody_system.nbody_to_si(
            system.mass.sum(), system.position[1].length()
        )

        instance = self.new_instance_of_an_optional_code(Tidymess, converter)
        assert instance is not None

        instance.parameters.tidal_model = 4
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 0
        instance.commit_parameters()

        self.assertAlmostRelativeEquals(
            instance.parameters.speed_of_light, constants.c
        )

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)

        instance.evolve_model(converter.to_si(1e3 | nbody_system.time))
        channel.copy()

        expected_position = [
            (9.2110173492271311e-2, 1.5422660672122340e-2, 0.0),
            (-2.3887569837471958e1, -3.9996654911555614, 0.0)
        ]

        for particle, (ex, ey, ez) in zip(instance.particles, expected_position):
            self.assertAlmostRelativeEquals(
                particle.x, converter.to_si(ex | nbody_system.length), places=4
            )
            self.assertAlmostRelativeEquals(
                particle.y, converter.to_si(ey | nbody_system.length), places=4
            )
            self.assertAlmostRelativeEquals(
                particle.z, converter.to_si(ez | nbody_system.length), places=4
            )

        expected_velocity = [
            (-4.5536563347452620e-4, 1.4689234497039288e-4, 0.0),
            (1.1809312651138797e-1, -3.8094610122321661e-2, 0.0)
        ]

        for particle, (evx, evy, evz) in zip(instance.particles, expected_velocity):
            self.assertAlmostRelativeEquals(
                particle.vx, converter.to_si(evx | nbody_system.speed), places=4
            )
            self.assertAlmostRelativeEquals(
                particle.vy, converter.to_si(evy | nbody_system.speed), places=4
            )
            self.assertAlmostRelativeEquals(
                particle.vz, converter.to_si(evz | nbody_system.speed), places=4
            )

        expected_spin = [
            (0.0, 0.0, 7.8439567544949240e-02),
            (0.0, 0.0, 3.8166572048686174e+00)
        ]

        for particle, (ewx, ewy, ewz) in zip(instance.particles, expected_spin):
            self.assertAlmostRelativeEquals(
                particle.wx, converter.to_si(ewx | 1/nbody_system.time), places=7
            )
            self.assertAlmostRelativeEquals(
                particle.wy, converter.to_si(ewy | 1/nbody_system.time), places=7
            )
            self.assertAlmostRelativeEquals(
                particle.wz, converter.to_si(ewz | 1/nbody_system.time), places=7
            )

        instance.stop()

    def test_nbody_units(self):
        """
        Evolve a figure 8 system in N-Body units with tides.

        For the expected values used to compare in this test,
        an identical simulation was ran in the Tidymess standalone
        package.
        """
        end_time = 2e3 | nbody_system.time
        system = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 4
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)

        instance.evolve_model(end_time)
        channel.copy()

        self.assertAlmostEquals(instance.model_time, end_time)

        expected_position = [
            (-8.3158513208432949e-1, -3.2772668917477010e-1, 0.0),
            (1.0473030710123414, 1.4115693325940498e-1, 0.0),
            (-2.1571793896223476e-1, 1.8656975594817923e-1, 0.0)
        ]

        for particle, (ex, ey, ez) in zip(instance.particles, expected_position):
            self.assertAlmostRelativeEquals(particle.x.number, ex, places=6)
            self.assertAlmostRelativeEquals(particle.y.number, ey, places=6)
            self.assertAlmostRelativeEquals(particle.z.number, ez, places=6)

        expected_velocity = [
            (-7.7344778284882154e-1, 2.9006761692492700e-1, 0.0),
            (-2.3236429722650778e-1, 4.6474040665908695e-1, 0.0),
            (1.0058120800752739, -7.5480802358400112e-1, 0.0)
        ]

        for particle, (evx, evy, evz) in zip(instance.particles, expected_velocity):
            self.assertAlmostRelativeEquals(particle.vx.number, evx, places=6)
            self.assertAlmostRelativeEquals(particle.vy.number, evy, places=6)
            self.assertAlmostRelativeEquals(particle.vz.number, evz, places=6)

        expected_spin = [
            (0.0, 0.0, -1.4197775319236833e-5),
            (0.0, 0.0, 2.1862179353662240e-7),
            (0.0, 0.0, -1.4059667978211637e-5)
        ]

        for particle, (ewx, ewy, ewz) in zip(instance.particles, expected_spin):
            self.assertAlmostRelativeEquals(particle.wx.number, ewx, places=5)
            self.assertAlmostRelativeEquals(particle.wy.number, ewy, places=5)
            self.assertAlmostRelativeEquals(particle.wz.number, ewz, places=5)

        # mass, radius, xi, kf, tau
        expected_attributes = [
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2)
        ]

        for particle, (m, r, xi, kf, tau) in zip(instance.particles, expected_attributes):
            self.assertAlmostEquals(particle.mass.number, m)
            self.assertAlmostEquals(particle.radius.number, r)
            self.assertAlmostEquals(particle.xi, xi)
            self.assertAlmostEquals(particle.kf, kf)
            self.assertAlmostEquals(particle.tau.number, tau)

        instance.stop()

    def test_tidal_model_1(self):
        """
        Evolve a system with tidal model 1.
        """
        end_time = 2e3 | nbody_system.time
        system = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 1
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)

        instance.evolve_model(end_time)
        channel.copy()

        self.assertAlmostEquals(instance.model_time, end_time)

        expected_position = [
            (-7.9995265821190509e-1, -3.3831974552195199e-1, 0.0),
            (1.0558216558327167, 1.2280654135281620e-1, 0.0),
            (-2.5586899764154514e-1, 2.1551320418599396e-1, 0.0)
        ]

        for particle, (ex, ey, ez) in zip(instance.particles, expected_position):
            self.assertAlmostRelativeEquals(particle.x.number, ex, places=6)
            self.assertAlmostRelativeEquals(particle.y.number, ey, places=6)
            self.assertAlmostRelativeEquals(particle.z.number, ez, places=6)

        expected_velocity = [
            (-8.3035491348610513e-1, 2.4451495581394561e-1, 0.0),
            (-1.9831541938121078e-1, 4.6603950585752918e-1, 0.0),
            (1.0286703328671680, -7.1055446167146374e-1, 0.0)
        ]

        for particle, (evx, evy, evz) in zip(instance.particles, expected_velocity):
            self.assertAlmostRelativeEquals(particle.vx.number, evx, places=6)
            self.assertAlmostRelativeEquals(particle.vy.number, evy, places=6)
            self.assertAlmostRelativeEquals(particle.vz.number, evz, places=6)

        expected_spin = [
            (0.0, 0.0, -5.6414636818972019e-16),
            (0.0, 0.0, 7.1678405859281773e-16),
            (0.0, 0.0, -2.8946120618110779e-16)
        ]

        for particle, (ewx, ewy, ewz) in zip(instance.particles, expected_spin):
            self.assertAlmostRelativeEquals(particle.wx.number, ewx, places=5)
            self.assertAlmostRelativeEquals(particle.wy.number, ewy, places=5)
            self.assertAlmostEquals(particle.wz.number, ewz)

        # mass, radius, xi, kf, tau
        expected_attributes = [
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2)
        ]

        for particle, (m, r, xi, kf, tau) in zip(instance.particles, expected_attributes):
            self.assertAlmostEquals(particle.mass.number, m)
            self.assertAlmostEquals(particle.radius.number, r)
            self.assertAlmostEquals(particle.xi, xi)
            self.assertAlmostEquals(particle.kf, kf)
            self.assertAlmostEquals(particle.tau.number, tau)

        instance.stop()

    def test_tidal_model_2(self):
        """
        Evolve a system with tidal model 2.
        """
        end_time = 2e3 | nbody_system.time
        system = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 2
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)

        instance.evolve_model(end_time)
        channel.copy()

        self.assertAlmostEquals(instance.model_time, end_time)

        expected_position = [
            (-8.3176224044049363e-1, -3.2766008362493232e-1, 0.0),
            (1.0472496444319142, 1.4126335057475434e-1, 0.0),
            (-2.1548740392829080e-1, 1.8639673306541499e-1, 0.0)
        ]

        for particle, (ex, ey, ez) in zip(instance.particles, expected_position):
            self.assertAlmostRelativeEquals(particle.x.number, ex, places=6)
            self.assertAlmostRelativeEquals(particle.y.number, ey, places=6)
            self.assertAlmostRelativeEquals(particle.z.number, ez, places=6)

        expected_velocity = [
            (-7.7311529372059096e-1, 2.9031053101895304e-1, 0.0),
            (-2.3256699425097455e-1, 4.6473080258147820e-1, 0.0),
            (1.0056822879716858, -7.5504133360041215e-1, 0.0)
        ]

        for particle, (evx, evy, evz) in zip(instance.particles, expected_velocity):
            self.assertAlmostRelativeEquals(particle.vx.number, evx, places=6)
            self.assertAlmostRelativeEquals(particle.vy.number, evy, places=6)
            self.assertAlmostRelativeEquals(particle.vz.number, evz, places=6)

        expected_spin = [
            (0.0, 0.0, -1.4343682408248869e-5),
            (0.0, 0.0, 2.2588894429400000e-7),
            (0.0, 0.0, -1.4232806418627201e-5)
        ]

        for particle, (ewx, ewy, ewz) in zip(instance.particles, expected_spin):
            self.assertAlmostRelativeEquals(particle.wx.number, ewx, places=5)
            self.assertAlmostRelativeEquals(particle.wy.number, ewy, places=5)
            self.assertAlmostEquals(particle.wz.number, ewz)

        # mass, radius, xi, kf, tau
        expected_attributes = [
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2)
        ]

        for particle, (m, r, xi, kf, tau) in zip(instance.particles, expected_attributes):
            self.assertAlmostEquals(particle.mass.number, m)
            self.assertAlmostEquals(particle.radius.number, r)
            self.assertAlmostEquals(particle.xi, xi)
            self.assertAlmostEquals(particle.kf, kf)
            self.assertAlmostEquals(particle.tau.number, tau)

        instance.stop()

    def test_tidal_model_3(self):
        """
        Evolve a system with tidal model 3.
        """
        end_time = 2e3 | nbody_system.time
        system = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 3
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)

        instance.evolve_model(end_time)
        channel.copy()

        self.assertAlmostEquals(instance.model_time, end_time)

        expected_position = [
            (-8.3158556071967982e-1, -3.2772652689780574e-1, 0.0),
            (1.0473029428899883, 1.4115718838003757e-1, 0.0),
            (-2.1571738128945700e-1, 1.8656933843546966e-1, 0.0)
        ]

        for particle, (ex, ey, ez) in zip(instance.particles, expected_position):
            self.assertAlmostRelativeEquals(particle.x.number, ex, places=6)
            self.assertAlmostRelativeEquals(particle.y.number, ey, places=6)
            self.assertAlmostRelativeEquals(particle.z.number, ez, places=6)

        expected_velocity = [
            (-7.7344697876688617e-1, 2.9006820552525181e-1, 0.0),
            (-2.3236478600092694e-1, 4.6474038388935018e-1, 0.0),
            (1.0058117647681448, -7.5480858941484674e-1, 0.0)
        ]

        for particle, (evx, evy, evz) in zip(instance.particles, expected_velocity):
            self.assertAlmostRelativeEquals(particle.vx.number, evx, places=6)
            self.assertAlmostRelativeEquals(particle.vy.number, evy, places=6)
            self.assertAlmostRelativeEquals(particle.vz.number, evz, places=6)

        expected_spin = [
            (0.0, 0.0, -1.4197809112719105e-5),
            (0.0, 0.0, 2.1862269039698131e-7),
            (0.0, 0.0, -1.4059672461331786e-5)
        ]

        for particle, (ewx, ewy, ewz) in zip(instance.particles, expected_spin):
            self.assertAlmostRelativeEquals(particle.wx.number, ewx, places=5)
            self.assertAlmostRelativeEquals(particle.wy.number, ewy, places=5)
            self.assertAlmostEquals(particle.wz.number, ewz)

        # mass, radius, xi, kf, tau
        expected_attributes = [
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2)
        ]

        for particle, (m, r, xi, kf, tau) in zip(instance.particles, expected_attributes):
            self.assertAlmostEquals(particle.mass.number, m)
            self.assertAlmostEquals(particle.radius.number, r)
            self.assertAlmostEquals(particle.xi, xi)
            self.assertAlmostEquals(particle.kf, kf)
            self.assertAlmostEquals(particle.tau.number, tau)

        instance.stop()

    def test_tidal_model_0(self):
        """
        Evolve a system with tidal model 0.
        """
        end_time = 2e3 | nbody_system.time
        system = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 0
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)

        instance.evolve_model(end_time)
        channel.copy()

        self.assertAlmostEquals(instance.model_time, end_time)

        expected_position = [
            (-7.9971116777386553e-1, -3.3839146769012957e-1, 0.0),
            (1.0558794343776374, 1.2267200440402883e-1, 0.0),
            (-2.5616826675029120e-1, 2.1571946325892974e-1, 0.0)
        ]

        for particle, (ex, ey, ez) in zip(instance.particles, expected_position):
            self.assertAlmostRelativeEquals(particle.x.number, ex, places=6)
            self.assertAlmostRelativeEquals(particle.y.number, ey, places=6)
            self.assertAlmostRelativeEquals(particle.z.number, ez, places=6)

        expected_velocity = [
            (-8.3077006303527656e-1, 2.4415103199669067e-1, 0.0),
            (-1.9807071425059053e-1, 4.6604655407076018e-1, 0.0),
            (1.0288407772857484, -7.1019758606749461e-1, 0.0)
        ]

        for particle, (evx, evy, evz) in zip(instance.particles, expected_velocity):
            self.assertAlmostRelativeEquals(particle.vx.number, evx, places=6)
            self.assertAlmostRelativeEquals(particle.vy.number, evy, places=6)
            self.assertAlmostRelativeEquals(particle.vz.number, evz, places=6)

        for particle in instance.particles:
            self.assertEquals(particle.wx.number, 0.0)
            self.assertEquals(particle.wy.number, 0.0)
            self.assertEquals(particle.wz.number, 0.0)

        # mass, radius, xi, kf, tau
        expected_attributes = [
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2),
            (1.0, 5e-2, 7e-2, 2e-2, 1e-2)
        ]

        for particle, (m, r, xi, kf, tau) in zip(instance.particles, expected_attributes):
            self.assertAlmostEquals(particle.mass.number, m)
            self.assertAlmostEquals(particle.radius.number, r)
            self.assertAlmostEquals(particle.xi, xi)
            self.assertAlmostEquals(particle.kf, kf)
            self.assertAlmostEquals(particle.tau.number, tau)

        instance.stop()

    def test_snapshot_dependency(self):
        """
        Test that varying diagnostic dt has no change on the result.
        This verifies that evolving to t_end in one `evolve_model`
        call gives the same answer as evolving to t_end incrementally.

        dt_diag controls how many times AMUSE creates a snapshot.
        """
        # System 1
        # --------
        t_end = 2e3 | nbody_system.time
        system1 = self.generate_figure8_system()

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 4
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system1)
        channel1 = instance.particles.new_channel_to(system1)

        instance.evolve_model(t_end)
        channel1.copy()

        instance.stop()

        # System 2
        # --------
        dt_diag = 1 | nbody_system.time
        system2 = self.generate_figure8_system()
        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.tidal_model = 4
        instance.parameters.dt_mode = 2
        instance.parameters.eta = 0.015625
        instance.parameters.initial_shape = 1
        instance.commit_parameters()

        instance.particles.add_particles(system2)
        channel2 = instance.particles.new_channel_to(system2)

        while instance.model_time < t_end:
            time = instance.model_time + dt_diag
            instance.evolve_model(time)
            channel2.copy()

        instance.stop()

        self.assertAlmostEquals(
            system1.position, system2.position
        )
        self.assertAlmostEquals(
            system1.velocity, system2.velocity, places=6
        )
        self.assertAlmostEquals(system1.wx, system2.wx)
        self.assertAlmostEquals(system1.wy, system2.wy)
        self.assertAlmostEquals(system1.wz, system2.wz)
        self.assertAlmostEquals(system1.mass, system2.mass)
        self.assertAlmostEquals(system1.radius, system2.radius)

    def test_stopping_conditions(self):
        """
        Test that collision detection works in Tidymess.

        A triple star system is initialized, with the expectation
        that stars 1 and 2 will collide.

        For the expected values, an identical standalone tidymess simulation
        was calculated.
        """
        p = Particles(3)
        p[0].name = 'Star 1'
        p[0].mass = 10 | nbody_system.mass
        p[0].radius = 1 | nbody_system.length
        p[0].x = -5 | nbody_system.length
        p[0].y = 0 | nbody_system.length
        p[0].z = 0 | nbody_system.length
        p[0].vx = 0 | nbody_system.speed
        p[0].vy = 0 | nbody_system.speed
        p[0].vz = 0 | nbody_system.speed
        p[0].xi = 0.0
        p[0].kf = 0.0
        p[0].tau = 0.0 | nbody_system.time
        p[0].wx = 0.0 | (1 / nbody_system.time)
        p[0].wy = 0.0 | (1 / nbody_system.time)
        p[0].wz = 0.0 | (1 / nbody_system.time)
        p[0].a_mb = 0 | nbody_system.time

        p[1].name = 'Star 2'
        p[1].mass = 10 | nbody_system.mass
        p[1].radius = 1 | nbody_system.length
        p[1].x = 5 | nbody_system.length
        p[1].y = 0 | nbody_system.length
        p[1].z = 0 | nbody_system.length
        p[1].vx = 0 | nbody_system.speed
        p[1].vy = 0 | nbody_system.speed
        p[1].vz = 0 | nbody_system.speed
        p[1].xi = 0.0
        p[1].kf = 0.0
        p[1].tau = 0.0 | nbody_system.time
        p[1].wx = 0.0 | (1 / nbody_system.time)
        p[1].wy = 0.0 | (1 / nbody_system.time)
        p[1].wz = 0.0 | (1 / nbody_system.time)
        p[1].a_mb = 0 | nbody_system.time

        p[2].name = 'Star 3'
        p[2].mass = 10 | nbody_system.mass
        p[2].radius = 1 | nbody_system.length
        p[2].x = -500 | nbody_system.length
        p[2].y = 0 | nbody_system.length
        p[2].z = 0 | nbody_system.length
        p[2].vx = 0 | nbody_system.speed
        p[2].vy = 0 | nbody_system.speed
        p[2].vz = 0 | nbody_system.speed
        p[2].xi = 0.0
        p[2].kf = 0.0
        p[2].tau = 0.0 | nbody_system.time
        p[2].wx = 0.0 | (1 / nbody_system.time)
        p[2].wy = 0.0 | (1 / nbody_system.time)
        p[2].wz = 0.0 | (1 / nbody_system.time)
        p[2].a_mb = 0 | nbody_system.time

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        cd = instance.stopping_conditions.collision_detection
        cd.enable()
        assert cd.is_supported() and cd.is_enabled()

        instance.parameters.collision_mode = 1
        instance.parameters.tidal_model = 0
        instance.parameters.dt_mode = 0
        instance.commit_parameters()

        instance.particles.add_particles(p)
        channel = instance.particles.new_channel_to(p)

        collision_hit = False

        dt_diag = 0.01 | nbody_system.time
        col_time = []
        while instance.model_time < 10 | nbody_system.time:
            time = instance.model_time + dt_diag
            instance.evolve_model(time)
            channel.copy()

            if cd.is_set():
                collision_hit = True
                col_time.append(time)
                break

        assert collision_hit
        col1 = cd.particles(0)
        col2 = cd.particles(1)

        assert col1[0] == p[0]
        assert col2[0] == p[1]

        assert len(col1) == 1
        assert len(col1) == len(col2)

        # check collision is defined correctly
        for pi, pj in zip(col1, col2):
            dr = pi.position - pj.position
            dist2 = (dr * dr).sum()
            rsum = pi.radius + pj.radius
            assert dist2 <= rsum**2

        # check that collision is expected from standalone tidymess
        expected_pos = [
            (-9.9264497812761021e-1, 0.0, 0.0),
            (9.9037034987767159e-1, 0.0, 0.0)
        ]
        for p, ep in zip((col1, col2), expected_pos):
            self.assertAlmostRelativeEquals(p.x[0].number, ep[0], places=7)
            self.assertAlmostRelativeEquals(p.y[0].number, ep[1], places=7)
            self.assertAlmostRelativeEquals(p.z[0].number, ep[2], places=7)

        expected_vel = [
            (2.0103752196648657, 0.0, 0.0),
            (-2.0109785401057181, 0.0, 0.0)
        ]
        for p, ev in zip((col1, col2), expected_vel):
            self.assertAlmostRelativeEquals(p.vx[0].number, ev[0], places=7)
            self.assertAlmostRelativeEquals(p.vy[0].number, ev[1], places=7)
            self.assertAlmostRelativeEquals(p.vz[0].number, ev[2], places=7)

        self.assertAlmostRelativeEquals(col_time[0].number, 7.5399999999998837, places=7)

        cd.disable()
        assert not cd.is_enabled()

        instance.stop()

    def test_triple_collision_stopping_conditions(self):
        """
        Initialize a triple star system on an equilateral triangle
        with vertices defined as:
            (0, 0, 0), (a, 0, 0), (a/2, a*sqrt(3)/2, 0)
        """
        a = 10 | nbody_system.length

        p = Particles(3)
        p[0].name = 'Star 1'
        p[0].mass = 10 | nbody_system.mass
        p[0].radius = 1 | nbody_system.length
        p[0].x = 0 | nbody_system.length
        p[0].y = 0 | nbody_system.length
        p[0].z = 0 | nbody_system.length
        p[0].vx = 0 | nbody_system.speed
        p[0].vy = 0 | nbody_system.speed
        p[0].vz = 0 | nbody_system.speed
        p[0].xi = 0.0
        p[0].kf = 0.0
        p[0].tau = 0.0 | nbody_system.time
        p[0].wx = 0.0 | (1 / nbody_system.time)
        p[0].wy = 0.0 | (1 / nbody_system.time)
        p[0].wz = 0.0 | (1 / nbody_system.time)
        p[0].a_mb = 0 | nbody_system.time

        p[1].name = 'Star 2'
        p[1].mass = 10 | nbody_system.mass
        p[1].radius = 1 | nbody_system.length
        p[1].x = a
        p[1].y = 0 | nbody_system.length
        p[1].z = 0 | nbody_system.length
        p[1].vx = 0 | nbody_system.speed
        p[1].vy = 0 | nbody_system.speed
        p[1].vz = 0 | nbody_system.speed
        p[1].xi = 0.0
        p[1].kf = 0.0
        p[1].tau = 0.0 | nbody_system.time
        p[1].wx = 0.0 | (1 / nbody_system.time)
        p[1].wy = 0.0 | (1 / nbody_system.time)
        p[1].wz = 0.0 | (1 / nbody_system.time)
        p[1].a_mb = 0 | nbody_system.time

        p[2].name = 'Star 3'
        p[2].mass = 10 | nbody_system.mass
        p[2].radius = 1 | nbody_system.length
        p[2].x = (a / 2)
        p[2].y = (a * (3**0.5) / 2)
        p[2].z = 0 | nbody_system.length
        p[2].vx = 0 | nbody_system.speed
        p[2].vy = 0 | nbody_system.speed
        p[2].vz = 0 | nbody_system.speed
        p[2].xi = 0.0
        p[2].kf = 0.0
        p[2].tau = 0.0 | nbody_system.time
        p[2].wx = 0.0 | (1 / nbody_system.time)
        p[2].wy = 0.0 | (1 / nbody_system.time)
        p[2].wz = 0.0 | (1 / nbody_system.time)
        p[2].a_mb = 0 | nbody_system.time

        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        cd = instance.stopping_conditions.collision_detection
        cd.enable()
        assert cd.is_supported() and cd.is_enabled()

        instance.parameters.collision_mode = 1
        instance.parameters.tidal_model = 0
        instance.parameters.dt_mode = 0
        instance.commit_parameters()

        instance.particles.add_particles(p)
        channel = instance.particles.new_channel_to(p)

        collision_hit = False

        dt_diag = 0.01 | nbody_system.time
        while instance.model_time < 10 | nbody_system.time:
            time = instance.model_time + dt_diag
            instance.evolve_model(time)
            channel.copy()

            if cd.is_set():
                collision_hit = True
                break

        assert collision_hit
        col1 = cd.particles(0)
        col2 = cd.particles(1)

        assert col1[0] == p[0]
        assert col1[1] == p[0]
        assert col1[2] == p[1]
        assert col2[0] == p[1]
        assert col2[1] == p[2]
        assert col2[2] == p[2]

        assert len(col1) == 3
        assert len(col2) == len(col1)

        # check collision is defined correctly
        for pi, pj in zip(col1, col2):
            dr = pi.position - pj.position
            dist2 = (dr * dr).sum()
            rsum = pi.radius + pj.radius
            assert dist2 <= rsum**2

        cd.disable()
        assert not cd.is_enabled()

        instance.stop()

    def test_warning_when_pn_set(self):
        """
        Test that a warning is generated if PN is set but not speed of light.
        """
        instance = self.new_instance_of_an_optional_code(Tidymess)
        assert instance is not None

        instance.parameters.pn_order = 1
        with pytest.warns(UserWarning):
            instance.commit_parameters()

        instance.stop()
