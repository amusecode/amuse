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

    def test_collision(self):
        print("Testing collision_detection")
        particles = datamodel.Particles(7)
        particles.mass = [1, 2.01, 4.02, 8.03, 16.04, 32.05, 64.06]* (0.000001 | nbody_system.mass)
        particles.radius = 0.01 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed
        
        instance = Petar(redirection="none")
        instance.parameters.dt_soft = 1/128 | nbody_system.time
        # instance.initialize_code()
        # instance.parameters.set_defaults()
        instance.parameters.epsilon_squared = 0 | nbody_system.length**2
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        
        for i in range(len(particles)):
            p = instance.particles[i]
            print(i, p.key, p.mass, p.x, p.vx)
        # print(instance.particles.x.value_in(nbody_system.length))
        for i in range(3):  # PH4 can handle only one collision (=closest) at a time
            for t in range(5):
                instance.evolve_model((t*0.2) | nbody_system.time)
                for i in range(len(particles)):
                    p = instance.particles[i]
                    print(i, p.key, p.mass, p.x, p.vx)
            # print(instance.particles.x.value_in(nbody_system.length))
            
            self.assertTrue(collisions.is_set())
            self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
            self.assertEqual(len(collisions.particles(0)), 1)
            self.assertEqual(len(collisions.particles(1)), 1)
            self.assertEqual(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 5 - i)
            self.assertEqual(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius), True)
        
            sticky_merged = datamodel.Particles(len(collisions.particles(0)))
            sticky_merged.mass = collisions.particles(0).mass + collisions.particles(1).mass
            sticky_merged.radius = collisions.particles(0).radius
            for p1, p2, merged in zip(collisions.particles(0), collisions.particles(1), sticky_merged):
                merged.position = (p1 + p2).center_of_mass()
                merged.velocity = (p1 + p2).center_of_mass_velocity()
        
            print(instance.model_time)
            print(instance.particles)
            instance.particles.remove_particles(collisions.particles(0) + collisions.particles(1))
            instance.particles.add_particles(sticky_merged)
        
        instance.evolve_model(1.0 | nbody_system.time)
        print()
        print(instance.model_time)
        print(instance.particles)
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 1.0 | nbody_system.time)
        self.assertEqual(len(collisions.particles(0)), 1)
        self.assertEqual(len(collisions.particles(1)), 1)
        self.assertEqual(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 2)
        self.assertEqual(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()
