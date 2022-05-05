from amuse.test.amusetest import TestWithMPI
import numpy

from amuse.community.petar.interface import PetarInterface, Petar

from amuse.units import nbody_system
from amuse import datamodel
from amuse.ic import plummer
from amuse.ic.plummer import new_plummer_model

from .gd_tests import _TestGravitationalDynamicsInterface


class TestPetarInterface(_TestGravitationalDynamicsInterface, TestWithMPI):
    PRECISION = 6

    def gravity_code_interface(self):
        return PetarInterface

    def reference_includes(self):
        return "Wang"

    def starting_particle_index(self):
        return 1

    def check_for_energy(self, *args, **kwargs):
        return self.assertAlmostEqual(*args, **kwargs, places=self.PRECISION)

    def test_reversed_time_allowed(self):
        self.skip("not supported")


class TestPetar(TestWithMPI):
    def test_small_plummer_model(self):
        particles = plummer.new_plummer_model(31)

        instance = Petar(number_of_workers=1)  # , debugger="xterm")
        instance.initialize_code()
        instance.particles.add_particles(particles)

        instance.evolve_model(0.1 | nbody_system.time)
        instance.synchronize_model()
        expected_positions = instance.particles.position
        instance.stop()
        positions_per_workers = []
        for n in [2, 3, 4, 5]:
            instance = Petar(number_of_workers=n)
            instance.initialize_code()
            instance.particles.add_particles(particles)

            instance.evolve_model(0.1 | nbody_system.time)
            instance.synchronize_model()
            positions_per_workers.append(instance.particles.position)
            instance.stop()

        for index, n in enumerate([2, 3, 4, 5]):
            self.assertAlmostEqual(
                expected_positions, positions_per_workers[index], 15)

    def test_reset_code(self):
        particles = new_plummer_model(50)
        particles.scale_to_standard()
        instance = Petar()
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 50)

        instance.reset()

        self.assertEqual(len(instance.particles), 0)
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 50)
        instance.stop()

    def test_adding_particles(self):
        petar = Petar()

        particles = datamodel.Particles(10)
        particles.position = ([0, 0, 0]) | nbody_system.length
        particles.velocity = ([1, 0, 0]) | nbody_system.speed
        particles.radius = 0 | nbody_system.length
        particles.mass = 0.1 | nbody_system.mass
        particles.x = numpy.linspace(1, 10, 10) | nbody_system.length
        particles.vx = numpy.linspace(1, 5, 10) | nbody_system.speed

        petar.particles.add_particles(particles)

        request = petar.particles.get_values_in_store_async(None, ["x"])
        request.wait()
        print(request.result())
        self.assertEqual(request.result()[0], particles.x)
        request = petar.particles.get_values_in_store_async(None, [
                                                              "x", "vx"])
        request.wait()
        print(request.result())
        self.assertEqual(request.result()[0], particles.x)
        self.assertEqual(request.result()[1], particles.vx)
        p = particles.copy()
        channel = petar.particles.new_channel_to(p)
        p.x = 0 | nbody_system.length
        p.vx = 0 | nbody_system.speed
        request = channel.copy_attributes_async(("x", "vx",), async_get=True)
        request.wait()
        self.assertEqual(p.x, particles.x)
        self.assertEqual(p.vx, particles.vx)
        p.x = 0 | nbody_system.length
        p.vx = 0 | nbody_system.speed
        channel = p.new_channel_to(petar.particles)
        request = channel.copy_attributes_async(
            ("x", "y", "z", "vx", "vy", "vz"), async_get=False, async_set=True)
        request.wait()
        self.assertEqual(p.x, petar.particles.x)
        self.assertEqual(p.vx, petar.particles.vx)
        channel = p.new_channel_to(particles)
        request = channel.copy_attributes_async(
            ("x", "y", "z", "vx", "vy", "vz"), async_get=False, async_set=True)
        request.wait()
        self.assertEqual(p.x, particles.x)
        self.assertEqual(p.vx, particles.vx)
        request = channel.copy_attributes_async(
            ("x", "y", "z", "vx", "vy", "vz"), async_get=True, async_set=False)
        request.wait()
        self.assertEqual(p.x, particles.x)
        self.assertEqual(p.vx, particles.vx)
