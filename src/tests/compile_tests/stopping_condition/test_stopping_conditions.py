from pathlib import Path

from amuse.support.testing.amusetest import TestWithMPI
from compile_tests.stopping_condition.interface import (
        ForTestingInterface, ForTestingInterfaceFortranModule)

# cello
from amuse.community.interface.stopping_conditions import StoppingConditions
from amuse.support.interface import InCodeComponentImplementation

from amuse.units import units
from amuse import datamodel
from amuse.support.exceptions import AmuseException
# from amuse.rfi.core import *
from amuse.community import NO_UNIT


class ForTesting(InCodeComponentImplementation):
    def __init__(self, exefile, **options):
        if 'community_interface' in options:
            interface = options['community_interface']
        else:
            interface = ForTestingInterface
        self.stopping_conditions = StoppingConditions(self)
        InCodeComponentImplementation.__init__(self, interface(exefile, **options), **options)
        self.my_particles = datamodel.Particles()

    def define_methods(self, object):
        self.stopping_conditions.define_methods(object)

    def new_particle(self, mass):
        particles = datamodel.Particles(len(mass))
        particles.mass = mass
        self.my_particles.add_particles(particles)
        return list(range(len(self.my_particles)-len(mass), len(self.my_particles)))

    def get_mass(self, indices):
        return self.my_particles.mass[indices]

    def delete_particle(self, particle):
        self.my_particles.remove_particle(particle)

    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_getter('particles', 'get_mass', names=("mass",))
        self.stopping_conditions.define_particle_set(object)


class TestInterface(TestWithMPI):

    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'c_worker')

    def test1(self):
        # ~ print self.exefile
        instance = ForTestingInterface(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        next = instance.next_index_for_stopping_condition()
        instance.stop()
        self.assertEqual(next, 1)

    def test2(self):
        instance = ForTesting(self.exefile)  # , debugger = "xterm")
        instance.initialize_code()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_supported())
        self.assertTrue(instance.stopping_conditions.collision_detection.is_supported())
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        instance.stop()

    def test3(self):
        instance = ForTesting(self.exefile)  # , debugger = "xterm")
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.enable()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.disable()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stop()

    def test4(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        # ~ print next,instance.stopping_conditions.pair_detection.type
        instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)

        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.stop()

    def test5(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())

        instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, 11)
        instance.set_stopping_condition_particle_index(next, 1, 12)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEqual(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEqual(12, instance.get_stopping_condition_particle_index(next, 1))
        instance.stop()

    def test6(self):
        instance = ForTesting(self.exefile)
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        instance.set_stopping_condition_info(next, instance.stopping_conditions.out_of_box_detection.type)
        self.assertTrue(instance.stopping_conditions.out_of_box_detection.is_set())
        instance.stop()

    def test7(self):
        instance = ForTesting(self.exefile)
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        instance.reset_stopping_conditions()

        pairs = [(11, 12), (0, 4), (3, 18), (7, 2)]
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, pairs[0][0])
        instance.set_stopping_condition_particle_index(next, 1, pairs[0][1])
        self.assertEqual(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEqual(12, instance.get_stopping_condition_particle_index(next, 1))
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(0)), 1)
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(1)), 1)

        for index1, index2 in pairs[1:]:
            next = instance.next_index_for_stopping_condition()
            instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)
            instance.set_stopping_condition_particle_index(next, 0, index1)
            instance.set_stopping_condition_particle_index(next, 1, index2)
            self.assertEqual(index1, instance.get_stopping_condition_particle_index(next, 0))
            self.assertEqual(index2, instance.get_stopping_condition_particle_index(next, 1))
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(0)), 4)
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(1)), 4)
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(2)), 0)

        self.assertEqual(instance.stopping_conditions.pair_detection.particles(0).mass,
            [first + 1 for first, second in pairs] | units.kg)
        self.assertEqual(instance.stopping_conditions.pair_detection.particles(1).mass,
            [second + 1 for first, second in pairs] | units.kg)
        instance.stop()

    def test8(self):
        instance = ForTesting(self.exefile)
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        self.assertRaises(AmuseException, instance.stopping_conditions.escaper_detection.enable, expected_message="Can't enable stopping condition 'escaper_detection', since 'ForTesting' does not support this condition.")
        instance.stop()

    def test9(self):
        instance = ForTestingInterface(self.exefile)
        instance.reset_stopping_conditions()
        nmax = 2048
        for i in range(nmax):
            next = instance.next_index_for_stopping_condition()
            # ~ print i, next
            self.assertEqual(next, i)
        instance.stop()


class TestInterfaceMP(TestWithMPI):

    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'c_mpi_ext_worker')

    def get_number_of_workers(self):
        return 3

    def test1(self):
        number_of_workers = 4
        instance = ForTestingInterface(self.exefile, number_of_workers=number_of_workers)
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()
        instance.enable_stopping_condition(1)
        nmax = 50
        for i in range(nmax):
            next = instance.next_index_for_stopping_condition()
            self.assertEqual(next, i)
        i, error = instance.get_number_of_stopping_conditions_set()
        self.assertEqual(error, 0)
        self.assertEqual(i, nmax)
        instance.mpi_collect_stopping_conditions()
        i, error = instance.get_number_of_stopping_conditions_set()
        self.assertEqual(error, 0)
        self.assertEqual(i, number_of_workers * nmax)

        instance.stop()

    def test2(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_distribute_stopping_conditions()

        # ~ print pair_detection.type
        instance.fire_condition(
            pair_detection.type,
            1, 2, -1
        )
        instance.mpi_collect_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEqual(len(pair_detection.particles(0)), self.get_number_of_workers())
        self.assertEqual(len(pair_detection.particles(1)), self.get_number_of_workers())
        self.assertEqual(pair_detection.particles(0).key, particles[1].key)
        self.assertEqual(pair_detection.particles(1).key, particles[2].key)
        self.assertEqual(pair_detection.particles(0).mass, [2, 2, 2] | units.kg)
        self.assertEqual(pair_detection.particles(1).mass, [3, 3, 3] | units.kg)
        instance.stop()

    def test5(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_distribute_stopping_conditions()
        for rank in range(self.get_number_of_workers()):
            # ~ print pair_detection.type
            instance.fire_condition(
                pair_detection.type,
                1, 2, rank
            )
            instance.mpi_collect_stopping_conditions()
            self.assertTrue(pair_detection.is_set())
            self.assertEqual(len(pair_detection.particles(0)), 1)
            self.assertEqual(len(pair_detection.particles(1)), 1)
            self.assertEqual(pair_detection.particles(0).key, particles[1].key)
            self.assertEqual(pair_detection.particles(1).key, particles[2].key)
            self.assertEqual(pair_detection.particles(0).mass, [2] | units.kg)
            self.assertEqual(pair_detection.particles(1).mass, [3] | units.kg)
            instance.reset_stopping_conditions()
            instance.stopping_conditions.pair_detection.enable()

        instance.stop()

    def test3(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_distribute_stopping_conditions()

        instance.fire_condition(
            pair_detection.type,
            1, 2, 0
        )
        instance.fire_condition(
            pair_detection.type,
            3, 4, 1
        )
        instance.fire_condition(
            pair_detection.type,
            5, 6, 2
        )
        instance.mpi_collect_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEqual(len(pair_detection.particles(0)), 3)
        self.assertEqual(len(pair_detection.particles(1)), 3)
        self.assertEqual(pair_detection.particles(0).key[0], particles[1].key)
        self.assertEqual(pair_detection.particles(1).key[0], particles[2].key)
        self.assertEqual(pair_detection.particles(0).key[1], particles[3].key)
        self.assertEqual(pair_detection.particles(1).key[1], particles[4].key)
        self.assertEqual(pair_detection.particles(0).key[2], particles[5].key)
        self.assertEqual(pair_detection.particles(1).key[2], particles[6].key)
        instance.reset_stopping_conditions()
        instance.stopping_conditions.pair_detection.enable()

        instance.stop()

    def test4(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_collect_stopping_conditions()

        instance.fire_condition(
            pair_detection.type,
            -1, -1, -1
        )
        instance.mpi_distribute_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEqual(len(pair_detection.particles(0)), 0)

        instance.stop()


class _TestInterfaceFortranSingleProcess:
    def get_number_of_workers(self):
        return 1

    def test1(self):
        instance = ForTestingInterface(self.exefile, number_of_workers=self.get_number_of_workers())
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        next = instance.next_index_for_stopping_condition()
        instance.stop()
        self.assertEqual(next, 1)

    def test2(self):
        instance = ForTesting(self.exefile, number_of_workers=self.get_number_of_workers())  # , debugger = "xterm")
        instance.initialize_code()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_supported())
        self.assertTrue(instance.stopping_conditions.collision_detection.is_supported())
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        instance.stop()

    def test3(self):
        instance = ForTesting(self.exefile, number_of_workers=self.get_number_of_workers())  # , debugger = "xterm")
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.enable()
        self.assertTrue(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stopping_conditions.pair_detection.disable()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_enabled())
        instance.stop()

    def test4(self):
        instance = ForTesting(self.exefile, number_of_workers=self.get_number_of_workers())
        instance.reset_stopping_conditions()

        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())

        instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)

        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.stop()

    def test5(self):
        instance = ForTesting(self.exefile, number_of_workers=self.get_number_of_workers())
        instance.reset_stopping_conditions()
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())

        instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, 11)
        instance.set_stopping_condition_particle_index(next, 1, 12)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEqual(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEqual(12, instance.get_stopping_condition_particle_index(next, 1))
        instance.stop()

    def test6(self):
        instance = ForTesting(self.exefile, number_of_workers=self.get_number_of_workers())
        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)
        instance.reset_stopping_conditions()

        pairs = [(11, 12), (0, 4), (3, 18), (7, 2)]
        next = instance.next_index_for_stopping_condition()
        self.assertFalse(instance.stopping_conditions.pair_detection.is_set())
        instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)
        instance.set_stopping_condition_particle_index(next, 0, pairs[0][0])
        instance.set_stopping_condition_particle_index(next, 1, pairs[0][1])
        self.assertEqual(11, instance.get_stopping_condition_particle_index(next, 0))
        self.assertEqual(12, instance.get_stopping_condition_particle_index(next, 1))
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(0)), 1)
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(1)), 1)

        for index1, index2 in pairs[1:]:
            next = instance.next_index_for_stopping_condition()
            instance.set_stopping_condition_info(next, instance.stopping_conditions.pair_detection.type)
            instance.set_stopping_condition_particle_index(next, 0, index1)
            instance.set_stopping_condition_particle_index(next, 1, index2)
            self.assertEqual(index1, instance.get_stopping_condition_particle_index(next, 0))
            self.assertEqual(index2, instance.get_stopping_condition_particle_index(next, 1))
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(0)), 4)
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(1)), 4)
        self.assertEqual(len(instance.stopping_conditions.pair_detection.particles(2)), 0)

        self.assertEqual(instance.stopping_conditions.pair_detection.particles(0).mass,
            [first + 1 for first, second in pairs] | units.kg)
        self.assertEqual(instance.stopping_conditions.pair_detection.particles(1).mass,
            [second + 1 for first, second in pairs] | units.kg)
        instance.stop()

    def test8(self):
        instance = ForTesting(self.exefile, number_of_workers=self.get_number_of_workers())
        instance.initialize_code()
        self.assertFalse(instance.stopping_conditions.escaper_detection.is_supported())
        self.assertRaises(AmuseException, instance.stopping_conditions.escaper_detection.enable, expected_message="Can't enable stopping condition 'escaper_detection', since 'ForTesting' does not support this condition.")
        instance.stop()

    def test9(self):
        instance = ForTestingInterface(self.exefile, number_of_workers=self.get_number_of_workers())
        instance.initialize_code()
        instance.reset_stopping_conditions()
        nmax = 2048
        for i in range(nmax):
            next = instance.next_index_for_stopping_condition()
            # ~ print i, next
            self.assertEqual(next, i)
        instance.stop()


class TestInterfaceFortran(TestWithMPI, _TestInterfaceFortranSingleProcess):
    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'f_worker')


class TestInterfaceFortranModule(TestWithMPI, _TestInterfaceFortranSingleProcess):
    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'f_ext_worker')


class TestInterfaceFortranModuleMultiprocess(TestWithMPI):
    @classmethod
    def setup_class(cls):
        cls.exefile = str(Path(__file__).parent / 'f_mpi_ext_worker')

    def get_number_of_workers(self):
        return 3

    def test1(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_distribute_stopping_conditions()

        # ~ print pair_detection.type
        instance.fire_condition(
            pair_detection.type,
            1, 2, -1
        )
        instance.mpi_collect_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEqual(len(pair_detection.particles(0)), self.get_number_of_workers())
        self.assertEqual(len(pair_detection.particles(1)), self.get_number_of_workers())
        self.assertEqual(pair_detection.particles(0).key, particles[1].key)
        self.assertEqual(pair_detection.particles(1).key, particles[2].key)
        self.assertEqual(pair_detection.particles(0).mass, [2, 2, 2] | units.kg)
        self.assertEqual(pair_detection.particles(1).mass, [3, 3, 3] | units.kg)
        instance.stop()

    def test2(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_distribute_stopping_conditions()
        for rank in range(self.get_number_of_workers()):
            # ~ print pair_detection.type
            instance.fire_condition(
                pair_detection.type,
                1, 2, rank
            )
            instance.mpi_collect_stopping_conditions()
            self.assertTrue(pair_detection.is_set())
            self.assertEqual(len(pair_detection.particles(0)), 1)
            self.assertEqual(len(pair_detection.particles(1)), 1)
            self.assertEqual(pair_detection.particles(0).key, particles[1].key)
            self.assertEqual(pair_detection.particles(1).key, particles[2].key)
            self.assertEqual(pair_detection.particles(0).mass, [2] | units.kg)
            self.assertEqual(pair_detection.particles(1).mass, [3] | units.kg)
            instance.reset_stopping_conditions()
            instance.stopping_conditions.pair_detection.enable()

        instance.stop()

    def test3(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_distribute_stopping_conditions()

        instance.fire_condition(
            pair_detection.type,
            1, 2, 0
        )
        instance.fire_condition(
            pair_detection.type,
            3, 4, 1
        )
        instance.fire_condition(
            pair_detection.type,
            5, 6, 2
        )
        instance.mpi_collect_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEqual(len(pair_detection.particles(0)), 3)
        self.assertEqual(len(pair_detection.particles(1)), 3)
        self.assertEqual(pair_detection.particles(0).key[0], particles[1].key)
        self.assertEqual(pair_detection.particles(1).key[0], particles[2].key)
        self.assertEqual(pair_detection.particles(0).key[1], particles[3].key)
        self.assertEqual(pair_detection.particles(1).key[1], particles[4].key)
        self.assertEqual(pair_detection.particles(0).key[2], particles[5].key)
        self.assertEqual(pair_detection.particles(1).key[2], particles[6].key)
        instance.reset_stopping_conditions()
        instance.stopping_conditions.pair_detection.enable()

        instance.stop()

    def test4(self):
        instance = ForTesting(
            self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=self.get_number_of_workers()
        )
        instance.initialize_code()
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()

        pair_detection = instance.stopping_conditions.pair_detection

        particles = datamodel.Particles(20)
        particles.mass = range(1, 21) | units.kg
        instance.particles.add_particles(particles)

        instance.stopping_conditions.pair_detection.enable()

        instance.mpi_collect_stopping_conditions()

        instance.fire_condition(
            pair_detection.type,
            -1, -1, -1
        )
        instance.mpi_distribute_stopping_conditions()
        self.assertTrue(pair_detection.is_set())
        self.assertEqual(len(pair_detection.particles(0)), 0)

        instance.stop()

    def test5(self):
        number_of_workers = 4
        instance = ForTestingInterface(self.exefile,
            community_interface=ForTestingInterfaceFortranModule,
            number_of_workers=number_of_workers)
        instance.reset_stopping_conditions()
        instance.mpi_setup_stopping_conditions()
        instance.enable_stopping_condition(1)
        nmax = 50
        for i in range(nmax):
            next = instance.next_index_for_stopping_condition()
            self.assertEqual(next, i)
        i, error = instance.get_number_of_stopping_conditions_set()
        self.assertEqual(error, 0)
        self.assertEqual(i, nmax)
        instance.mpi_collect_stopping_conditions()
        i, error = instance.get_number_of_stopping_conditions_set()
        self.assertEqual(error, 0)
        self.assertEqual(i, number_of_workers * nmax)

        instance.stop()
