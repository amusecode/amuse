from equilibrium_tides_averaged import EqTides, make_tsunami_ic
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import numpy as np
import pytest
import tsunami

from amuse.datamodel import Particle, Particles
from amuse.ext.orbital_elements import generate_binaries, orbital_elements
from amuse.support.testing.amusetest import TestWithMPI
from amuse.community.tsunami.interface import Tsunami, TsunamiInterface
from amuse.units import nbody_system as ns, units as u

# if True, show comparison plots for certain tests
show_plots = False

class TestTsunamiInterface(TestWithMPI):

    def test_initialization(self):
        """Test Tsunami initialization."""
        instance = self.new_instance_of_an_optional_code(TsunamiInterface)
        assert instance is not None

        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test_setters_and_getters(self):
        """Test TsunamiInterface setters and getters."""
        instance = self.new_instance_of_an_optional_code(TsunamiInterface)
        assert instance is not None

        self.assertEquals(0, instance.initialize_code())

        instance.set_units(2, 2)
        self.assertEquals(instance.get_Mscale()['Mscale'], 2)
        self.assertEquals(instance.get_Lscale()['Lscale'], 2)
        self.assertAlmostRelativeEquals(
            instance.get_Vscale()['Vscale'], 29.78469182967693
        )
        self.assertAlmostRelativeEquals(
            instance.get_Tscale()['Tscale'], 0.31831589802181726
        )
        self.assertAlmostRelativeEquals(
            instance.get_speed_of_light()['speed_of_light'],
            10065.320121972596
        )

        instance.set_equilibrium_tides(True)
        self.assertEquals(
            instance.get_equilibrium_tides()['equilibrium_tides'],
            True
        )

        instance.set_dynamical_tides(True)
        self.assertEquals(
            instance.get_dynamical_tides()['dynamical_tides'],
            True
        )

        instance.set_alpha(2)
        self.assertEquals(instance.get_alpha()['alpha'], 2)
        instance.set_beta(2)
        self.assertEquals(instance.get_beta()['beta'], 2)
        instance.set_gamma(2)
        self.assertEquals(instance.get_gamma()['gamma'], 2)

        instance.set_pn(True)
        self.assertEquals(instance.get_pn()['pn'], True)
        instance.set_pn1(True)
        self.assertEquals(instance.get_pn1()['pn1'], True)
        instance.set_pn2(True)
        self.assertEquals(instance.get_pn2()['pn2'], True)
        instance.set_pn25(True)
        self.assertEquals(instance.get_pn25()['pn25'], True)
        instance.set_pn3(True)
        self.assertEquals(instance.get_pn3()['pn3'], True)
        instance.set_pn35(True)
        self.assertEquals(instance.get_pn35()['pn35'], True)

        instance.set_units(1, 1)

        self.assertEquals(0, instance.commit_parameters())

        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 0)

        result = instance.new_particle(3,0,1,3,0,0,0,0,0,0,0)
        self.assertEquals(result['index_of_the_particle'], 0)

        result = instance.new_particle(4,0,-2,-1,0,0,0,0,1,1,1)
        self.assertEquals(result['index_of_the_particle'], 1)

        result = instance.new_particle(5,0,1,-1,0,0,0,0,0,0,0)
        self.assertEquals(result['index_of_the_particle'], 2)

        instance.commit_particles()
        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 3)

        result = instance.get_state(0)
        self.assertEquals(result['mass'], 3.0)
        self.assertEquals(result['radius'], 0.0)
        self.assertEquals(result['x'], 1.0)
        self.assertEquals(result['y'], 3.0)
        self.assertEquals(result['z'], 0.0)
        self.assertEquals(result['vx'], 0.0)
        self.assertEquals(result['vy'], 0.0)
        self.assertEquals(result['vz'], 0.0)
        self.assertEquals(result['wx'], 0.0)
        self.assertEquals(result['wy'], 0.0)
        self.assertEquals(result['wz'], 0.0)

        result = instance.get_state(1)
        self.assertEquals(result['mass'], 4)
        self.assertEquals(result['radius'], 0)
        self.assertEquals(result['x'], -2)
        self.assertEquals(result['y'], -1)
        self.assertEquals(result['z'], 0)
        self.assertEquals(result['vx'], 0)
        self.assertEquals(result['vy'], 0)
        self.assertEquals(result['vz'], 0)
        self.assertEquals(result['wx'], 1)
        self.assertEquals(result['wy'], 1)
        self.assertEquals(result['wz'], 1)

        result = instance.get_state(2)
        self.assertEquals(result['mass'], 5)
        self.assertEquals(result['radius'], 0)
        self.assertEquals(result['x'], 1)
        self.assertEquals(result['y'], -1)
        self.assertEquals(result['z'], 0)
        self.assertEquals(result['vx'], 0)
        self.assertEquals(result['vy'], 0)
        self.assertEquals(result['vz'], 0)
        self.assertEquals(result['wx'], 0)
        self.assertEquals(result['wy'], 0)
        self.assertEquals(result['wz'], 0)

        instance.set_mass(0, 1.5)
        self.assertEquals(instance.get_mass(0)['mass'], 1.5)

        self.assertEquals(instance.get_radius(0)['radius'], 0.0)

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

        result = instance.get_spin(1)
        self.assertEquals(result['wx'], 1)
        self.assertEquals(result['wy'], 1)
        self.assertEquals(result['wz'], 1)

        instance.delete_particle(0)
        instance.recommit_particles()
        result = instance.get_number_of_particles()
        self.assertEquals(result['number_of_particles'], 2)

        result = instance.get_state(1)
        self.assertEquals(result['mass'], 4)
        self.assertEquals(result['radius'], 0)
        self.assertEquals(result['x'], -2)
        self.assertEquals(result['y'], -1)
        self.assertEquals(result['z'], 0)
        self.assertEquals(result['vx'], 0)
        self.assertEquals(result['vy'], 0)
        self.assertEquals(result['vz'], 0)
        self.assertEquals(result['wx'], 1)
        self.assertEquals(result['wy'], 1)
        self.assertEquals(result['wz'], 1)

        result = instance.get_state(2)
        self.assertEquals(result['mass'], 5)
        self.assertEquals(result['radius'], 0)
        self.assertEquals(result['x'], 1)
        self.assertEquals(result['y'], -1)
        self.assertEquals(result['z'], 0)
        self.assertEquals(result['vx'], 0)
        self.assertEquals(result['vy'], 0)
        self.assertEquals(result['vz'], 0)
        self.assertEquals(result['wx'], 0)
        self.assertEquals(result['wy'], 0)
        self.assertEquals(result['wz'], 0)

        instance.set_state(
            1,
            5, 5,
            5, 5, 5,
            5, 5, 5,
            5, 5, 5
        )
        result = instance.get_state(1)
        self.assertEquals(result['mass'], 5)
        self.assertEquals(result['x'], 5)
        self.assertEquals(result['y'], 5)
        self.assertEquals(result['z'], 5)
        self.assertEquals(result['vx'], 5)
        self.assertEquals(result['vy'], 5)
        self.assertEquals(result['vz'], 5)

        # check that this is still the same
        result = instance.get_state(2)
        self.assertEquals(result['mass'], 5)
        self.assertEquals(result['radius'], 0)
        self.assertEquals(result['x'], 1)
        self.assertEquals(result['y'], -1)
        self.assertEquals(result['z'], 0)
        self.assertEquals(result['vx'], 0)
        self.assertEquals(result['vy'], 0)
        self.assertEquals(result['vz'], 0)
        self.assertEquals(result['wx'], 0)
        self.assertEquals(result['wy'], 0)
        self.assertEquals(result['wz'], 0)

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()


class TestTsunami(TestWithMPI):
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
        HD80606[0].wx = 0.0 | 1 / u.s
        HD80606[0].wy = 0.0 | 1 / u.s
        HD80606[0].wz = 2.97188607137e-6 | 1 / u.s

        HD80606[1].name = 'HD80606b'
        HD80606[1].mass = planet.mass
        HD80606[1].radius = 68488.3446 | u.km
        HD80606[1].x = planet.x.as_quantity_in(u.km)
        HD80606[1].y = planet.y.as_quantity_in(u.km)
        HD80606[1].z = planet.z.as_quantity_in(u.km)
        HD80606[1].vx = planet.vx.as_quantity_in(u.kms)
        HD80606[1].vy = planet.vy.as_quantity_in(u.kms)
        HD80606[1].vz = planet.vz.as_quantity_in(u.kms)
        HD80606[1].wx = 0.0 | 1 / u.s
        HD80606[1].wy = 0.0 | 1 / u.s
        HD80606[1].wz = 0.000145444104333 | 1 / u.s

        return HD80606

    def generate_pythagorean(self):
        """
        Generate 3 particles in a pythagorean triangle configuration.
        Initial conditions courtesy of Dr. Alessandro Alberto Trani.
        """
        p = Particles(3)

        p[0].mass = 3 | ns.mass
        p[0].radius = 0 | ns.length
        p[0].x = 1 | ns.length
        p[0].y = 3 | ns.length
        p[0].z = 0 | ns.length
        p[0].vx = 0 | ns.speed
        p[0].vy = 0 | ns.speed
        p[0].vz = 0 | ns.speed
        p[0].wx = 0 | 1 / ns.time
        p[0].wy = 0 | 1 / ns.time
        p[0].wz = 0 | 1 / ns.time

        p[1].mass = 4 | ns.mass
        p[1].radius = 0 | ns.length
        p[1].x = -2 | ns.length
        p[1].y = -1 | ns.length
        p[1].z = 0 | ns.length
        p[1].vx = 0 | ns.speed
        p[1].vy = 0 | ns.speed
        p[1].vz = 0 | ns.speed
        p[1].wx = 0 | 1 / ns.time
        p[1].wy = 0 | 1 / ns.time
        p[1].wz = 0 | 1 / ns.time

        p[2].mass = 5 | ns.mass
        p[2].radius = 0 | ns.length
        p[2].x = 1 | ns.length
        p[2].y = -1 | ns.length
        p[2].z = 0 | ns.length
        p[2].vx = 0 | ns.speed
        p[2].vy = 0 | ns.speed
        p[2].vz = 0 | ns.speed
        p[2].wx = 0 | 1 / ns.time
        p[2].wy = 0 | 1 / ns.time
        p[2].wz = 0 | 1 / ns.time

        return p

    def generate_earth_moon_initial_conditions(self):
        """
        Earth and moon initial conditions taken from JPL Ephemeris.

        JPL Ephemeris, Pasadena, USA, Horizons System

            * Start time        : A.D. 2026-Jan-01 00:00:00.0000 TDB
            * Stop  time        : A.D. 2026-Jan-01 00:01:00.0000 TDB
            * Coordinate Center : Solar System Barycenter
        """
        p = Particles(2)

        p[0].name = "Earth"
        p[0].mass = 5.97219e24 | u.kg
        p[0].radius = 6371.01 | u.km
        p[0].xi = 0.3308
        p[0].kf = 0.933
        p[0].tau = 60. | u.s
        p[0].wx = 0.0 | u.rad/u.s
        p[0].wy = 0.0 | u.rad/u.s
        p[0].wz = 0.0 | u.rad/u.s
        p[0].x = -2.653100241556548E+07 | u.km
        p[0].y = 1.439468995740296E+08 | u.km
        p[0].z = 1.080681311843544E+04 | u.km
        p[0].vx = -2.977650610770464E+01 | u.kms
        p[0].vy = -5.395962660572101E+00 | u.kms
        p[0].vz = 1.753836198843395E-04 | u.kms

        p[1].name = "Moon"
        p[1].mass = 7.349e22 | u.kg
        p[1].radius = 1737.53 | u.km
        p[1].xi = 0.394
        p[1].kf = 0.4
        p[1].tau = 60. | u.s
        p[1].wx = 0.0 | u.rad/u.s
        p[1].wy = 0.0 | u.rad/u.s
        p[1].wz = 0.0 | u.rad/u.s
        p[1].x = -2.638667668807356E+07 | u.km
        p[1].y = 1.442762954049252E+08 | u.km
        p[1].z = 4.255978946162760E+04 | u.km
        p[1].vx = -3.078082024537686E+01 | u.kms
        p[1].vy = -4.975097451815738E+00 | u.kms
        p[1].vz = 5.760594917691098E-03 | u.kms

        p.move_to_center()

        return p

    def _validate_tsunami_state_equality(self, state: list, particle: Particle) -> None:
        """
        Validate that a state retrieved by `Tsunami.get_state`
        matches the state of `particle` exactly.

        Parameters
        ----------
        state : list
            List containing the state of the particle,
            as returned by `get_state`.
        particle : amuse.datamodel.particles.particle
            Particle to validate against `state`.

        Raises
        ------
        AssertionError :
            If the states don't match.
        """
        self.assertEquals(state[0], particle.mass)
        self.assertEquals(state[1], particle.radius)
        self.assertEquals(state[2], particle.x)
        self.assertEquals(state[3], particle.y)
        self.assertEquals(state[4], particle.z)
        self.assertEquals(state[5], particle.vx)
        self.assertEquals(state[6], particle.vy)
        self.assertEquals(state[7], particle.vz)
        self.assertEquals(state[8], particle.wx)
        self.assertEquals(state[9], particle.wy)
        self.assertEquals(state[10], particle.wz)


    def _validate_tsunami_state_relative_equality(
        self,
        state: list,
        particle: Particle,
        places: int = 7
    ) -> None:
        """
        Validate that a state retrieved via `Tsunami.get_state` matches
        the corresponding attributes of `particle`, within relative tolerance.

        Parameters
        ----------
        state : list
            State vector as returned by `get_state`, ordered as
            [mass, radius, x, y, z, vx, vy, vz, wx, wy, wz].
        particle : amuse.datamodel.particles.Particle
            Particle instance whose attributes are compared against `state`.
        places : int, optional, default=7
            Number of decimal places of relative precision required for equality.

        Raises
        ------
        AssertionError :
            If any corresponding value in `state` and `particle` does not match
            within `places` relative precision.
        """
        self.assertAlmostRelativeEquals(state[0], particle.mass, places=places)
        self.assertAlmostRelativeEquals(state[1], particle.radius, places=places)
        self.assertAlmostRelativeEquals(state[2], particle.x, places=places)
        self.assertAlmostRelativeEquals(state[3], particle.y, places=places)
        self.assertAlmostRelativeEquals(state[4], particle.z, places=places)
        self.assertAlmostRelativeEquals(state[5], particle.vx, places=places)
        self.assertAlmostRelativeEquals(state[6], particle.vy, places=places)
        self.assertAlmostRelativeEquals(state[7], particle.vz, places=places)
        self.assertAlmostRelativeEquals(state[8], particle.wx, places=places)
        self.assertAlmostRelativeEquals(state[9], particle.wy, places=places)
        self.assertAlmostRelativeEquals(state[10], particle.wz, places=places)

    def test_add_and_delete_particles(self):
        """Test adding and deleting particles in Tsunami."""
        system = self.generate_pythagorean()

        instance = self.new_instance_of_an_optional_code(Tsunami)
        assert instance is not None

        instance.particles.add_particles(system)
        instance.commit_particles()

        self.assertEquals(instance.get_number_of_particles(), len(system))

        state0 = instance.get_state(0)
        state1 = instance.get_state(1)
        state2 = instance.get_state(2)

        self._validate_tsunami_state_equality(state0, system[0])
        self._validate_tsunami_state_equality(state1, system[1])
        self._validate_tsunami_state_equality(state2, system[2])

        self.assertEquals(system[0].key, instance.particles[0].key)
        self.assertEquals(system[1].key, instance.particles[1].key)
        self.assertEquals(system[2].key, instance.particles[2].key)

        instance.particles.remove_particle(system[0])
        self.assertEquals(instance.get_number_of_particles(), len(system)-1)

        state1 = instance.get_state(1)
        state2 = instance.get_state(2)

        self._validate_tsunami_state_equality(state1, system[1])
        self._validate_tsunami_state_equality(state2, system[2])

        instance.stop()

    def test_tsunami_pythagorean_triple(self):
        """
        Test identical to `test_tsunami.py` from the TSUNAMI package.
        Simulation results should be identical to standalone TSUNAMI.
        """
        system = self.generate_pythagorean()

        instance = self.new_instance_of_an_optional_code(Tsunami, redirection='none')
        assert instance is not None

        instance.parameters.wPNs = False
        instance.parameters.wEqTides = False
        instance.parameters.wDynTides = False
        instance.commit_parameters()

        instance.particles.add_particles(system)
        instance.commit_particles()

        t_end = 65 | ns.time
        dt = 0.1 | ns.time
        while instance.model_time < t_end:
            instance.evolve_model(instance.model_time + dt)

        # validate energies
        self.assertAlmostRelativeEquals(instance.potential_energy.number, 22.461771188652676)
        self.assertAlmostRelativeEquals(instance.kinetic_energy.number, 9.645104521979146)
        self.assertAlmostRelativeEquals(instance.deltaE, 3.055333763768898e-13)

        self.assertAlmostRelativeEquals(instance.model_time, t_end, places=2)
        state0 = instance.get_state(0)
        state1 = instance.get_state(1)
        state2 = instance.get_state(2)

        # compare with standalone TSUNAMI simulation results
        self.assertEquals(state0[0], 3 | ns.mass)
        self.assertEquals(state0[1], 0 | ns.length)
        self.assertAlmostRelativeEquals(state0[2], 4.15498855 | ns.length, places=7)
        self.assertAlmostRelativeEquals(state0[3], 12.01338805 | ns.length, places=7)
        self.assertEquals(state0[4], 0 | ns.length)
        self.assertAlmostRelativeEquals(state0[5], 0.59219296 | ns.speed, places=7)
        self.assertAlmostRelativeEquals(state0[6], 1.75695995 | ns.speed, places=7)
        self.assertEquals(state0[7], 0 | ns.speed)
        self.assertEquals(state0[8], 0 | 1 / ns.time)
        self.assertEquals(state0[9], 0 | 1 / ns.time)
        self.assertEquals(state0[10], 0 | 1 / ns.time)

        self.assertEquals(state1[0], 4 | ns.mass)
        self.assertEquals(state1[1], 0 | ns.length)
        self.assertAlmostRelativeEquals(state1[2], -0.86025151 | ns.length, places=7)
        self.assertAlmostRelativeEquals(state1[3], -4.09450317 | ns.length, places=7)
        self.assertEquals(state1[4], 0 | ns.length)
        self.assertAlmostRelativeEquals(state1[5], 0.61291976 | ns.speed, places=7)
        self.assertAlmostRelativeEquals(state1[6], -0.92151518 | ns.speed, places=7)
        self.assertEquals(state1[7], 0 | ns.speed)
        self.assertEquals(state1[8], 0 | 1 / ns.time)
        self.assertEquals(state1[9], 0 | 1 / ns.time)
        self.assertEquals(state1[10], 0 | 1 / ns.time)

        self.assertEquals(state2[0], 5 | ns.mass)
        self.assertEquals(state2[1], 0 | ns.length)
        self.assertAlmostRelativeEquals(state2[2], -1.80479192 | ns.length, places=7)
        self.assertAlmostRelativeEquals(state2[3], -3.9324303 | ns.length, places=7)
        self.assertEquals(state2[4], 0 | ns.length)
        self.assertAlmostRelativeEquals(state2[5], -0.84565158 | ns.speed, places=7)
        self.assertAlmostRelativeEquals(state2[6], -0.31696382 | ns.speed, places=7)
        self.assertEquals(state2[7], 0 | ns.speed)
        self.assertEquals(state2[8], 0 | 1 / ns.time)
        self.assertEquals(state2[9], 0 | 1 / ns.time)
        self.assertEquals(state2[10], 0 | 1 / ns.time)

        instance.stop()

    def test_delete_particle_and_restart(self):
        """
        This test simulates a pythagorean triangle system
        identical to `test_tsunami_pythagorean_triangle`,
        but instead of evolving to time t=65, it first evolves
        to time t=30, then deletes a particle, then adds it back
        in before finishing the integration to t=65.
        Both this test and `test_tsunami_pythagorean_triangle`
        should agree. Due to integration differences, the results
        of the two different simulations differ by more than 1 decimal
        point, and the states are not compared in this test.
        Set `show_plots=True` at the top of the test file to show the
        comparison plot.
        """
        system = self.generate_pythagorean()

        instance = self.new_instance_of_an_optional_code(Tsunami, redirection='none')
        assert instance is not None

        instance.parameters.wPNs = False
        instance.parameters.wEqTides = False
        instance.parameters.wDynTides = False
        instance.commit_parameters()

        instance.particles.add_particles(system)
        channel = instance.particles.new_channel_to(system)
        instance.commit_particles()

        particles = []

        t_midpoint = 30 | ns.time
        dt = 0.1 | ns.time
        while instance.model_time < t_midpoint:
            instance.evolve_model(instance.model_time + dt)
            channel.copy()
            particles.append(system.copy())

        self.assertAlmostRelativeEquals(instance.model_time, t_midpoint, places=2)

        # save original state after evolving to t=30
        pars = instance.particles.copy()

        # remove particle 0
        instance.particles.remove_particle(instance.particles[0])
        self.assertEquals(instance.get_number_of_particles(), 2)

        # add 'particle 0' back in
        instance.particles.add_particle(pars[0])
        instance.recommit_particles()
        self.assertEquals(instance.get_number_of_particles(), 3)

        # check that new starting state is the same as t=30 state before deleting
        # particle 0 now has id=3 because ids are unique and increasing monotonically
        state0 = instance.get_state(3)
        state1 = instance.get_state(1)
        state2 = instance.get_state(2)
        self._validate_tsunami_state_relative_equality(state0, pars[0], places=10)
        self._validate_tsunami_state_relative_equality(state1, pars[1], places=10)
        self._validate_tsunami_state_relative_equality(state2, pars[2], places=10)

        t_end = 65 | ns.time
        dt = 0.1 | ns.time
        while instance.model_time < t_end:
            instance.evolve_model(instance.model_time + dt)
            channel.copy()
            particles.append(system.copy())

        self.assertAlmostRelativeEquals(instance.model_time, t_end, places=2)

        if show_plots:
            colors = ['#483d8b', '#d81b60', '#dbb0ff']
            ax = self.plot_particles_xy(
                particles, colors=colors, labels=['star1', 'star2', 'star3']
            )
            self.plot_tsunami_pythagorean_triangle(ax)
            ax.legend(edgecolor='w')
            plt.show()

        instance.stop()

    def test_tsunami_does_not_rescale_to_com_frame(self):
        """
        By default Tsunami will rescale a system to its center
        of mass frame. Test that this behavior is turned off.
        """
        p = Particles(3)
        p.mass = [1,4,8] | ns.mass
        p.radius = [1,1,2] | ns.length
        p.position = [[0,1,2],[3,4,5],[6,7,8]] | ns.length
        p.velocity = [[6,7,8],[9,10,11],[12,13,14]] | ns.speed
        p.wx = [12,15,19] | 1 / ns.time
        p.wy = [13,16,20] | 1 / ns.time
        p.wz = [14,17,21] | 1 / ns.time

        instance = self.new_instance_of_an_optional_code(Tsunami, redirection='none')
        instance.particles.add_particles(p)
        instance.commit_particles()

        state0 = instance.get_state(0)
        state1 = instance.get_state(1)
        state2 = instance.get_state(2)

        self._validate_tsunami_state_equality(state0, p[0])
        self._validate_tsunami_state_equality(state1, p[1])
        self._validate_tsunami_state_equality(state2, p[2])

        instance.delete_particle(1)

        state0 = instance.get_state(0)
        state2 = instance.get_state(2)
        self._validate_tsunami_state_relative_equality(state0, p[0])
        self._validate_tsunami_state_relative_equality(state2, p[2])

        instance.set_mass(0, 4 | ns.mass)
        instance.set_position(0, 67 | ns.length, 68 | ns.length, 69 | ns.length)
        instance.set_velocity(2, 67 | ns.speed, 68 | ns.speed, 69 | ns.speed)

        self.assertEquals(instance.get_mass(0), 4 | ns.mass)
        self.assertEquals(
            instance.get_position(0),
            [67 | ns.length, 68 | ns.length, 69 | ns.length]
        )
        self.assertEquals(
            instance.get_velocity(2),
            [67 | ns.speed, 68 | ns.speed, 69 | ns.speed]
        )
        p3 = Particle(
            mass=1 | ns.mass,
            radius=100 | ns.length,
            position=[89,90,91] | ns.length,
            velocity=[92,93,94] | ns.speed,
            wx=95 | 1 / ns.time,
            wy=96 | 1 / ns.time,
            wz=97 | 1 / ns.time
        )
        instance.particles.add_particle(p3)

        state3 = instance.get_state(3)
        self._validate_tsunami_state_equality(state3, p3)

        # check that particles 0 and 2 have not changed
        state0 = instance.get_state(0)
        state2 = instance.get_state(2)
        self._validate_tsunami_state_relative_equality(state0, instance.particles[0])
        self._validate_tsunami_state_relative_equality(state2, instance.particles[2])

        instance.stop()

    def test_converter(self):
        """
        Test that the converter works to convert physical units
        to nbody units.

        NOTE : To recreate a Tsunami simulation, please setup
        the AMUSE converter as follows:
            `converter = nbody_system.nbody_to_si(M | u.MSun, L | u.au)`
        Where `M` and `L` match the values used in Tsunami,
        ie. `code = tsunami.Tsunami(M, L)`.

        Due to small unit / constant value differences (ie. G),
        do not expect AMUSE Tsunami and standalone Tsunami to
        match perfectly.
        """
        p = self.generate_earth_moon_initial_conditions()
        converter = ns.nbody_to_si(1 | u.MSun, 1 | u.au)

        instance = self.new_instance_of_an_optional_code(
            Tsunami, convert_nbody=converter
        )
        assert instance is not None
        instance.particles.add_particles(p)

        nbody_pos_amuse = converter.to_nbody(instance.particles.position)
        nbody_vel_amuse = converter.to_nbody(instance.particles.velocity)

        p2 = self.generate_earth_moon_initial_conditions()
        m = np.ascontiguousarray(p2.mass.value_in(u.Msun))
        r = np.ascontiguousarray(p2.radius.value_in(u.AU))
        pos = np.ascontiguousarray(p2.position.value_in(u.AU))
        vel = np.ascontiguousarray(p2.velocity.value_in(u.AU / u.yr))
        wx = np.ascontiguousarray(p2.wx.value_in(1 / u.yr))
        wy = np.ascontiguousarray(p2.wy.value_in(1 / u.yr))
        wz = np.ascontiguousarray(p2.wz.value_in(1 / u.yr))
        spin = np.vstack([wx,wy,wz]).T

        code = tsunami.Tsunami(1, 1)
        code.Conf.wExt = True

        # convert tsunami ics to nbody units
        m_nb = m / code.Mscale
        r_nb = r / code.Lscale
        pos_nb = pos / code.Lscale
        vel_nb = vel / (code.Lscale / code.Tscale)
        spin_nb = np.ascontiguousarray(spin * code.Tscale)
        pt = np.ones_like(m_nb, dtype=np.int64) * -1

        code.add_particle_set(pos_nb, vel_nb, m_nb, r_nb, pt, spin_nb)
        code.sync_internal_state(pos_nb, vel_nb, spin_nb)

        # check that the converter works
        self.assertAlmostRelativeEquals(nbody_pos_amuse.number, pos_nb)
        self.assertAlmostRelativeEquals(nbody_vel_amuse.number, vel_nb, places=3)

        instance.stop()

    def test_earth_moon_system_physical_units(self):
        """
        Test the Earth Moon system with physical units.
        Because of slight differences between the values
        of units / constants used in Tsunami and AMUSE,
        it is not possible to perfectly recreate Tsunami
        simulation in AMUSE when physical units are used.
        """
        # SETUP AMUSE RUN
        # ---------------
        p = self.generate_earth_moon_initial_conditions()
        converter = ns.nbody_to_si(1 | u.MSun, 1 | u.au)

        instance = self.new_instance_of_an_optional_code(
            Tsunami, convert_nbody=converter
        )
        assert instance is not None

        instance.particles.add_particles(p)

        nbody_pos_amuse = converter.to_nbody(instance.particles.position)
        nbody_vel_amuse = converter.to_nbody(instance.particles.velocity)

        end_time = (2 | u.yr).as_quantity_in(u.s)
        dt = 50000 | u.s

        # EVOLVE AMUSE
        # ------------
        while instance.model_time < end_time:
            instance.evolve_model(instance.model_time + dt)

        # SETUP IDENTICAL TSUNAMI RUN
        # ---------------------------
        p2 = self.generate_earth_moon_initial_conditions()
        m = np.ascontiguousarray(p2.mass.value_in(u.Msun))
        r = np.ascontiguousarray(p2.radius.value_in(u.AU))
        pos = np.ascontiguousarray(p2.position.value_in(u.AU))
        vel = np.ascontiguousarray(p2.velocity.value_in(u.AU / u.yr))
        wx = np.ascontiguousarray(p2.wx.value_in(1 / u.yr))
        wy = np.ascontiguousarray(p2.wy.value_in(1 / u.yr))
        wz = np.ascontiguousarray(p2.wz.value_in(1 / u.yr))
        spin = np.vstack([wx,wy,wz]).T

        code = tsunami.Tsunami(1, 1)
        code.Conf.wExt = True

        # convert tsunami ics to nbody units
        m_nb = m / code.Mscale
        r_nb = r / code.Lscale
        pos_nb = pos / code.Lscale
        vel_nb = vel / (code.Lscale / code.Tscale)
        spin_nb = np.ascontiguousarray(spin * code.Tscale)
        pt = np.ones_like(m_nb, dtype=np.int64) * -1

        code.add_particle_set(pos_nb, vel_nb, m_nb, r_nb, pt, spin_nb)
        code.sync_internal_state(pos_nb, vel_nb, spin_nb)

        # check that the converter works
        self.assertAlmostRelativeEquals(nbody_pos_amuse.number, pos_nb)
        self.assertAlmostRelativeEquals(nbody_vel_amuse.number, vel_nb, places=3)

        t_end_yr = 2
        dt_yr = dt.value_in(u.yr)
        t_end_nb = t_end_yr / code.Tscale
        dt_nb = dt_yr / code.Tscale
        t = 0

        # EVOLVE TSUNAMI
        # --------------
        while t < t_end_nb:
            code.evolve_system(t + dt_nb)
            t = code.time
            code.sync_internal_state(pos_nb, vel_nb, spin_nb)

        # VALIDATE THAT BOTH RUNS MATCH
        # -----------------------------
        pos_au = instance.particles.position.as_quantity_in(u.AU)
        vel_kms = instance.particles.velocity.as_quantity_in(u.kms)
        wx = instance.particles.wx.as_quantity_in(1/u.yr)
        wy = instance.particles.wy.as_quantity_in(1/u.yr)
        wz = instance.particles.wz.as_quantity_in(1/u.yr)

        tsunami_pos = pos_nb * code.Lscale | u.AU
        tsunami_vel = vel_nb * code.Vscale | u.kms
        tsunami_spin = spin_nb / code.Tscale | 1 / u.yr

        self.assertAlmostRelativeEquals(pos_au, tsunami_pos, places=2)

        self.assertAlmostRelativeEquals(vel_kms[0,0], tsunami_vel[0,0], places=2)
        self.assertAlmostRelativeEquals(vel_kms[0,1], tsunami_vel[0,1], places=2)
        self.assertAlmostRelativeEquals(vel_kms[1,0], tsunami_vel[1,0], places=2)
        self.assertAlmostRelativeEquals(vel_kms[1,1], tsunami_vel[1,1], places=2)

        # as of right now, Tsunami and AMUSE do not use the same exact unit values,
        # namely the value of G and/or the value of pi. thus it is not possible
        # to reproduce the exact values
        with pytest.raises(AssertionError):
            self.assertAlmostRelativeEquals(vel_kms[0,2], tsunami_vel[0,2], places=2)
            self.assertAlmostRelativeEquals(vel_kms[1,2], tsunami_vel[1,2], places=2)

        self.assertEquals(wx[0], tsunami_spin[0,0])
        self.assertEquals(wx[1], tsunami_spin[1,0])
        self.assertEquals(wy[0], tsunami_spin[0,1])
        self.assertEquals(wy[1], tsunami_spin[1,1])
        self.assertEquals(wz[0], tsunami_spin[0,2])
        self.assertEquals(wz[1], tsunami_spin[1,2])

        instance.stop()

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
        converter = ns.nbody_to_si(
            system.mass.sum(), 0.455 | u.AU
        )

        instance = self.new_instance_of_an_optional_code(
            Tsunami, convert_nbody=converter
        )
        assert instance is not None

        instance.particles.add_particles(system)
        instance.commit_particles()

        m = np.ascontiguousarray(system.mass.value_in(u.Msun))
        r = np.ascontiguousarray(system.radius.value_in(u.AU))
        pos = np.ascontiguousarray(system.position.value_in(u.AU))
        vel = np.ascontiguousarray(system.velocity.value_in(u.AU / u.yr))
        wx = np.ascontiguousarray(system.wx.value_in(1 / u.yr))
        wy = np.ascontiguousarray(system.wy.value_in(1 / u.yr))
        wz = np.ascontiguousarray(system.wz.value_in(1 / u.yr))
        spin = np.vstack([wx,wy,wz]).T
        code = tsunami.Tsunami(1,1)
        code.Conf.wExt = True
        m_nb = m
        r_nb = r
        pos_nb = pos
        vel_nb = vel / (2 * np.pi)
        spin_nb = np.ascontiguousarray(spin / (2 * np.pi))
        pt = np.ones_like(m_nb, dtype=np.int64) * -1
        code.add_particle_set(pos_nb, vel_nb, m_nb, r_nb, pt, spin_nb)
        code.sync_internal_state(pos_nb, vel_nb, spin_nb)

        positions = instance.particles.position.value_in(u.AU)
        velocities = instance.particles.velocity.value_in(u.km/u.s)
        wx = instance.particles.wx.value_in(u.rad/u.yr)
        wy = instance.particles.wy.value_in(u.rad/u.yr)
        wz = instance.particles.wz.value_in(u.rad/u.yr)

        self.assertAlmostRelativeEquals(positions, pos_nb * code.Lscale)
        self.assertAlmostRelativeEquals(velocities, vel_nb * code.Vscale, places=4)
        self.assertAlmostRelativeEquals(np.vstack([wx, wy, wz]).T, spin_nb / code.Tscale, places=4)

        instance.stop()

    def test_tides(self):
        """
        Test Tsunami tidal evolution. An AMUSE simulation
        is compared with a standalone Tsunami tidal evolution test
        script written by Dr. Alessandro Alberto Trani. The initial
        conditions are also from this test script.
        """
        ############################
        # Tsunami Initial Conditions
        ############################
        KU = tsunami.KeplerUtils()
        m1 = 1  # MSun
        m2 = 1  # MSun
        R1 = 1 * KU.RSun2au
        R2 = 1 * KU.RSun2au
        # apsidal motion constant
        kap = 0.05
        k1 = kap
        k2 = kap
        # time-lag, first in seconds, converted to N-body
        tausec = 1e3
        tau1 = tausec / KU.yr2sec / KU.Tscale
        tau2 = tausec / KU.yr2sec / KU.Tscale
        # gyration radius, given inertia I = M * (rg * R)^2
        rg = 0.28
        rg1 = rg
        rg2 = rg
        # initial orbit
        a0, e0 = 0.09, 0.3
        # spin, first in period (days) then in angular frequency
        Ps1 = 100 / 365.25 / KU.Tscale
        spin1 = 2 * np.pi / Ps1
        Ps2 = 0.0015 / KU.Tscale
        spin2 = 2 * np.pi / Ps2

        ET = EqTides(m1, m2, R1, R2, k1, k2, tau1, tau2, rg1, rg2, rotdist=False, usespin=True)
        ET.setup_orbital_vectors(a0, e0)
        obl1 = np.radians(30)
        obl2 = 0.0
        ET.setup_spin_vectors(spin1, spin2, obl1, obl2)
        ET.setup_vectors()
        ET.initialize_integrator()

        m, R, p, v, spin, _ = make_tsunami_ic(m1, m2, R1, R2, a0, e0, ET.spin1_vec, ET.spin2_vec)

        # final time, in units of initial orbital period (in N-body units)
        tfin = 5e5 * ET.P0
        dt = 5e1 * ET.P0

        kaps = np.array([k1, k2])
        polyind = np.array([0.0, 0.0])  # only for dynamical tide, null values
        gyrad = np.array([rg1, rg2])

        ##################
        # AMUSE SIMULATION
        ##################
        converter = ns.nbody_to_si(1 | u.MSun, 1 | u.au)
        instance = self.new_instance_of_an_optional_code(
            Tsunami, convert_nbody=converter
        )
        assert instance is not None
        # enable tides
        instance.parameters.wEqTides = True

        mass = m | u.MSun
        radii = R | u.au
        pos = p | u.au
        vel = converter.to_si(v | ns.speed).as_quantity_in(u.AU/u.yr)
        spin_si = converter.to_si(spin | 1 / ns.time).as_quantity_in(1/u.yr)
        tau = np.array([tausec, tausec]) | u.s

        p = Particles(2)
        p.mass = mass
        p.radius = radii
        p.position = pos
        p.velocity = vel
        p.spin = spin_si
        p[0].name = 'A'
        p[1].name = 'B'

        # add particles
        instance.particles.add_particles(p)
        channel = instance.particles.new_channel_to(p)
        instance.commit_particles()
        # initialize tides
        instance.initialize_tidal_parameters(
            kaps,
            tau,
            polyind,
            gyrad
        )

        # convert final time and dt to SI units
        tfin_si = converter.to_si(tfin | ns.time).as_quantity_in(u.yr)
        dt_si = converter.to_si(dt | ns.time).as_quantity_in(u.yr)

        pars = [p.copy()]
        times = [instance.model_time]
        while instance.model_time < tfin_si:
            instance.evolve_model(instance.model_time + dt_si)
            channel.copy()
            pars.append(p.copy())
            times.append(instance.model_time)

        elements = [orbital_elements(p) for p in pars]
        smas = [e[2].as_quantity_in(u.au) for e in elements]
        eccentricities = [e[3] for e in elements]
        print(eccentricities[-1])
        # check that last semi major axis & e measurements agrees with Tsunami
        self.assertAlmostRelativeEquals(
            smas[-1], 0.090719705156567 | u.au, places=5
        )
        self.assertAlmostRelativeEquals(
            eccentricities[-1], 0.2992015985318265400, places=5
        )
        # check that final times are the same
        self.assertAlmostRelativeEquals(converter.to_nbody(
            instance.model_time).number, 59983.559799160954, places=1
        )

        instance.stop()

    def plot_particles_xy(
        self,
        particles: list[Particles],
        colors: list[str],
        labels: list[str] | None = None,
        figsize: tuple[float, float] = (6,6)
    ) -> Axes:
        """
        Plot the xy position of a list of `Particles`.

        Parameters
        ----------
        particles : list[Particles]
            List of `amuse.datamodel.particles.Particles`, where each `Particles`
            in the list is a simulation snapshot containing each `Particle` in the system.
        colors : list[ColorType]
            List of colors for each particle. Should have same length as the number of particles
            in the system.
        labels : list[str]
            Label for each particle. Should have same length as the number of particles
            in the system.
        figsize : tuple[float, float]
            Size of the figure.

        Returns
        -------
        ax : Axes
        """
        data = [[], [], []]

        for p in particles:
            for i in range(len(p)):
                data[i].append((p[i].x.number, p[i].y.number))

        fig = plt.figure(figsize=figsize, tight_layout=True)
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.minorticks_on()
        ax.tick_params(
            axis='both', length=2, direction='in',
            which='both', right=True, top=True
        )

        for i, (color, pts) in enumerate(zip(colors, data)):
            xs, ys = zip(*pts)
            l = None if labels is None else labels[i]
            ax.plot(xs, ys, color=color, linewidth=0.8, label=l)

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        return ax

    def plot_tsunami_pythagorean_triangle(self, ax: Axes) -> None:
        """Plot `test_tsunami.py` from the standalone TSUNAMI package."""
        import tsunami
        code = tsunami.Tsunami(1.0, 1.0)

        code.Conf.wPNs = False
        code.Conf.wEqTides = False
        code.Conf.wDynTides = False
        code.Conf.wExt = False

        code.Conf.dcoll = 0.0

        m = np.array([3., 4., 5.])
        p = np.array([[1.,3.,0.], [-2.,-1.,0.], [1.,-1.,0.]])
        v = np.array([[0.,0.,0.], [0.,0.,0.], [0.,0.,0.]])
        R = np.array([0., 0., 0.])
        sp = np.array([[0.,0.,0.], [0.,0.,0.], [0.,0.,0.]])
        st = np.array([-1, -1, -1])
        code.add_particle_set(p, v, m, R, st, sp)

        code.sync_internal_state(p, v, sp)
        totp = [p.copy()]

        dt = 0.1
        ft = 65

        time = 0
        while time < ft:
            time = time + dt

            code.evolve_system(time)
            time = code.time
            code.sync_internal_state(p, v, sp)

            totp.append(p.copy())

        totp = np.vstack(totp)
        code.sync_masses(m)
        code.sync_radii(R)

        colors = ['#26dcba', '#7d7ff3', '#cfe23c']

        ax.plot(totp[::3,0], totp[::3,1], lw=0.8, ls='--', color=colors[0], label='star1 tsunami')
        ax.plot(totp[1::3,0], totp[1::3,1], lw=0.8, ls='--', color=colors[1], label='star2 tsunami')
        ax.plot(totp[2::3,0], totp[2::3,1], lw=0.8, ls='--', color=colors[2], label='star3 tsunami')
