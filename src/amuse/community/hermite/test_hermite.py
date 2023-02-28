from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.hermite.interface import HermiteInterface, Hermite

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.ic import plummer
from amuse.ic.plummer import new_plummer_model

from amuse.test.suite.codes_tests.gd_tests import (
    _TestGravitationalDynamicsInterface,
)
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestHermiteInterface(_TestGravitationalDynamicsInterface, TestWithMPI):
    def gravity_code_interface(self):
        return HermiteInterface

    def reference_includes(self):
        return "Hut"

    def test6(self):
        hermite = HermiteInterface()
        hermite.initialize_code()

        hermite.new_particle([10,10],[-1,1],[0,0], [0,0], [0,0], [0,0], [0,0], [1,1])
        retrieved_state = hermite.get_state(0)

        retr = hermite.get_potential_at_point(0.01, 0,0,0)
        self.assertEqual(retr['phi'], -20.0)
        hermite.cleanup_code()
        hermite.stop()

    def test7(self):
        instance = HermiteInterface()
        instance.initialize_code()
        instance.set_eps2(0.1 * 0.1)
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass=10.0, radius=1.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
        id2,errorcode = instance.new_particle(mass=10.0, radius=1.0, x=2.0, y=0.0, z=0.0, vx=10.0, vy=0.0, vz=0.0)

        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEqual(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2 + 0.1**2), 8)

        potential, errorcode = instance.get_potential(id2)
        self.assertEqual(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2 + 0.1**2), 8)

        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        instance.cleanup_code()
        instance.stop()

        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 10.0]) / 2.0)


    def test8(self):
        instance = HermiteInterface()
        instance.initialize_code()
        instance.set_eps2(0)
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass=10.0, radius=1.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
        id2,errorcode = instance.new_particle(mass=1.0, radius=1.0, x=2.0, y=0.0, z=0.0, vx=10.0, vy=0.0, vz=0.0)

        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEqual(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -1.0 / numpy.sqrt(2.0**2), 8)
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        instance.cleanup_code()
        instance.stop()

        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 1.0]) / 2.0)


    def test9(self):
        print("Test HermiteInterface evolve_model")
        instance = HermiteInterface()
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.set_dt_param(0.001))
        self.assertEqual(0, instance.set_end_time_accuracy_factor(0.0))
        self.assertEqual(0, instance.commit_parameters())

        # Set up an equal-mass binary on a circular orbit:
        self.assertEqual([0, 0], list(instance.new_particle(0.5,  0.5, 0, 0,  0, 0.5, 0, 0.01).values()))
        self.assertEqual([1, 0], list(instance.new_particle(0.5,  -0.5, 0, 0,  0,-0.5, 0, 0.01).values()))
        self.assertEqual(0, instance.commit_particles())

        self.assertEqual(0, instance.evolve_model(math.pi))
        for result, expected in zip(list(instance.get_position(0).values()), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(list(instance.get_position(1).values()), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)

        self.assertEqual(0, instance.evolve_model(2 * math.pi))
        for result, expected in zip(list(instance.get_position(0).values()), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(list(instance.get_position(1).values()), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)

        self.assertEqual(0, instance.cleanup_code())
        instance.cleanup_code()
        instance.stop()


class TestHermite(TestWithMPI):
    def new_system_of_sun_and_earth(self):
        stars = datamodel.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))

        return stars


    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = Hermite(convert_nbody)
        hermite.initialize_code()
        hermite.parameters.epsilon_squared = 0.0 | units.AU**2
        hermite.parameters.end_time_accuracy_factor = 0.0
        hermite.dt_dia = 5000

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        hermite.particles.add_particles(stars)

        hermite.evolve_model(365.0 | units.day)
        hermite.particles.copy_values_of_all_attributes_to(stars)

        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostRelativeEqual(position_at_start, position_after_full_rotation, 6)

        hermite.evolve_model(365.0 + (365.0 / 2) | units.day)

        hermite.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 3)

        hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)

        hermite.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 3)

        hermite.cleanup_code()
        hermite.stop()


    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Hermite(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.dt_dia = 5000

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.particles.add_particles(stars)

        for x in range(1, 500, 10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_all_attributes_to(stars)
            stars.savepoint()

        if HAS_MATPLOTLIB:
            figure = pyplot.figure()
            plot = figure.add_subplot(1,1,1)

            x_points = earth.get_timeline_of_attribute("x")
            y_points = earth.get_timeline_of_attribute("y")

            x_points_in_AU = [t_x[1].value_in(units.AU) for t_x in x_points]
            y_points_in_AU = [t_x1[1].value_in(units.AU) for t_x1 in y_points]

            plot.scatter(x_points_in_AU,y_points_in_AU, color="b", marker='o')

            plot.set_xlim(-1.5, 1.5)
            plot.set_ylim(-1.5, 1.5)


            test_results_path = self.get_path_to_results()
            output_file = os.path.join(test_results_path, "hermite-earth-sun2.svg")
            figure.savefig(output_file)



        instance.cleanup_code()
        instance.stop()

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Hermite(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | units.AU**2
        instance.dt_dia = 5000

        stars = datamodel.Stars(2)
        star1 = stars[0]
        star2 = stars[1]

        star1.mass = units.MSun(1.0)
        star1.position = units.AU(numpy.array((-1.0,0.0,0.0)))
        star1.velocity = units.AUd(numpy.array((0.0,0.0,0.0)))
        star1.radius = units.RSun(1.0)

        star2.mass = units.MSun(1.0)
        star2.position = units.AU(numpy.array((1.0,0.0,0.0)))
        star2.velocity = units.AUd(numpy.array((0.0,0.0,0.0)))
        star2.radius = units.RSun(100.0)

        instance.particles.add_particles(stars)

        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_all_attributes_to(stars)
            stars.savepoint()


        instance.cleanup_code()
        instance.stop()

    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Hermite(convert_nbody)
        instance.initialize_code()

        particles = datamodel.Particles(2)
        self.assertEqual(len(instance.particles), 0)

        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s


        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 2)

        instance.particles.mass =  [17.0, 33.0] | units.kg


        self.assertEqual(instance.get_mass(0), 17.0| units.kg) 
        self.assertEqual(instance.get_mass(1), 33.0| units.kg)  

        instance.stop()

    def test5(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Hermite(convert_nbody)
        instance.initialize_code()

        particles = datamodel.Particles(2)
        self.assertEqual(len(instance.particles), 0)

        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s


        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 2)
        instance.set_state(1, 16|units.kg, 
                           20.0|units.m, 40.0|units.m, 60.0|units.m, 
                           1.0|units.ms, 1.0|units.ms, 1.0|units.ms , 
                           20.0|units.m)

        curr_state =  instance.get_state(1)
        for expected, actual in zip([16|units.kg, 
                           20.0|units.m, 40.0|units.m, 60.0|units.m, 
                           1.0|units.ms, 1.0|units.ms, 1.0|units.ms , 
                           20.0|units.m], curr_state):
            self.assertAlmostRelativeEquals(expected, actual)
        instance.stop()

        self.assertEqual(curr_state[0], 16|units.kg, 8)

    def test6(self):
        print("Test6: Testing Hermite parameters")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = Hermite(convert_nbody)

        value = instance.get_eps2()
        self.assertEqual(0.0 | units.AU**2 , value)
        self.assertAlmostEqual(0.0 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEqual(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)

        value = instance.get_dt_param()
        self.assertEqual(0.03 , value)
        self.assertAlmostEqual(0.03 , instance.parameters.dt_param)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.dt_param = x
            self.assertAlmostEqual(x, instance.parameters.dt_param)

        value = instance.get_dt_dia()
        self.assertAlmostEqual(1.0 | units.yr, value)
        self.assertAlmostEqual(1.0 | units.yr, instance.parameters.dt_dia, in_units=units.yr)
        for x in [0.1, 10.0, 100.0]:
            instance.parameters.dt_dia = x | units.yr
            self.assertAlmostEqual(x | units.yr, instance.parameters.dt_dia, in_units=units.yr)

        value = instance.get_begin_time()
        self.assertEqual(0.0| units.yr, value)
        self.assertAlmostEqual(0.0 | units.yr, instance.parameters.begin_time, in_units=units.yr)
        for x in [1.0, 10.0, 100.0]:
            instance.parameters.begin_time = x | units.yr
            self.assertAlmostEqual(x | units.yr, instance.parameters.begin_time, in_units=units.yr)
        instance.stop()

    def test7(self):
        print("Test7: Testing effect of Hermite parameter epsilon_squared")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)

        particles = datamodel.Particles(2)
        sun = particles[0]
        sun.mass = 1.0 | units.MSun
        sun.position = [0.0, 0.0, 0.0] | units.AU
        sun.velocity = [0.0, 0.0, 0.0] | units.AU / units.yr
        sun.radius = 1.0 | units.RSun

        earth = particles[1]
        earth.mass = 5.9736e24 | units.kg
        earth.radius = 6371.0 | units.km
        earth.position = [0.0, 1.0, 0.0] | units.AU
        earth.velocity = [2.0*numpy.pi, -0.0001, 0.0] | units.AU / units.yr

        initial_direction = math.atan((earth.velocity[0]/earth.velocity[1]))
        final_direction = []
        for log_eps2 in range(-9,10,2):
            instance = Hermite(convert_nbody)
            instance.parameters.end_time_accuracy_factor = 0.0
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            final_direction.append(math.atan((instance.particles[1].velocity[0]/
                instance.particles[1].velocity[1])))
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEqual(abs(final_direction[0]), abs(initial_direction+math.pi/2.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEqual(final_direction[-1], initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(final_direction[i+1]-final_direction[i]) for i in range(len(final_direction)-1)]
        self.assertEqual(delta[len(final_direction)//2 -1], max(delta))

    def test8(self):
        print("Testing Hermite collision_detection")
        particles = datamodel.Particles(7)
        particles.mass = 0.001 | nbody_system.mass
        particles.radius = 0.01 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed

        instance = Hermite()
        instance.initialize_code()
        instance.parameters.set_defaults()
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        instance.evolve_model(1.0 | nbody_system.time)

        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
        self.assertEqual(len(collisions.particles(0)), 3)
        self.assertEqual(len(collisions.particles(1)), 3)
        self.assertEqual(len(particles - collisions.particles(0) - collisions.particles(1)), 1)
        self.assertEqual(abs(collisions.particles(0).x - collisions.particles(1).x) <= 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True, True, True])

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
        self.assertEqual(abs(collisions.particles(0).x - collisions.particles(1).x) <= 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()


    def test10(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Hermite(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.stopping_conditions_number_of_steps = 10
        self.assertEqual(instance.parameters.stopping_conditions_number_of_steps,10)

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        instance.particles.add_particles(stars)
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(365.0 | units.day)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        instance.particles.copy_values_of_all_attributes_to(stars)

        instance.cleanup_code()

        instance.stop()

    def test11(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        instance = Hermite()
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 2
        self.assertEqual(instance.parameters.stopping_conditions_number_of_steps, 2)
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody_system.time)

        instance.stop()

    def test12(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,1.00] | nbody_system.length
        particles.y = [0.0,0.0] | nbody_system.length
        particles.z = [0.0,0.0] | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  [5.1,0.0] | nbody_system.speed
        particles.vy =  [0.0,0.0] | nbody_system.speed
        particles.vz =  [0.0,0.0]| nbody_system.speed
        particles.mass = [0.1,0.1] | nbody_system.mass

        instance = Hermite()
        instance.initialize_code()
        instance.parameters.stopping_conditions_out_of_box_size = .5 | nbody_system.length
        self.assertEqual(instance.parameters.stopping_conditions_out_of_box_size, .5 | nbody_system.length)
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.out_of_box_detection.enable()
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.out_of_box_detection.is_set())
        self.assertAlmostEqual(instance.stopping_conditions.out_of_box_detection.particles(0).x, 1.0 |nbody_system.length, 3)
        instance.stop()

    def test13(self):
        particles = plummer.new_plummer_model(31)

        instance = Hermite(number_of_workers=1)#, debugger="xterm")
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.01 | nbody_system.length ** 2
        instance.particles.add_particles(particles)

        instance.evolve_model(0.1 | nbody_system.time)
        instance.synchronize_model()
        expected_positions = instance.particles.position
        instance.stop()
        positions_per_workers = []
        for n in [2,3,4,5]:
            instance = Hermite(number_of_workers=n)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 0.01 | nbody_system.length ** 2
            instance.particles.add_particles(particles)

            instance.evolve_model(0.1 | nbody_system.time)
            instance.synchronize_model()
            positions_per_workers.append(instance.particles.position)
            instance.stop()


        for index, n in enumerate([2,3,4,5]):
            self.assertAlmostEqual(expected_positions, positions_per_workers[index], 15)

    def test14(self):
        instance = Hermite()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | nbody_system.length**2

        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)

        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, 1.0 | nbody_system.length, zero, zero)
        self.assertAlmostEqual(fx, 0.0 | nbody_system.acceleration, 6)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration, 6)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration, 6)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            potential1 = instance.get_potential_at_point(zero, x1, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            fx1, fy1, fz1 = instance.get_gravity_at_point(zero, x1, zero, zero)

            self.assertAlmostEqual(fy0, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fz0, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fy1, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fz1, 0.0 | nbody_system.acceleration, 6)

            self.assertAlmostEqual(fx0, -1.0 * fx1, 5)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0, 2)
            self.assertAlmostEqual(potential0, potential1, 5)
        instance.stop()

    def test15(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0.0 | nbody_system.length
        particles.z = 0.0 | nbody_system.length
        particles.vx =  0.0 | nbody_system.speed
        particles.vy =  0.0 | nbody_system.speed
        particles.vz =  0.0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        instance = Hermite()
        instance.particles.add_particles(particles) 
        instance.commit_particles()
        self.assertEqual(instance.particles[0].radius, 0.0 | nbody_system.length)
        p = datamodel.Particle(
            x=1.0  | nbody_system.length,
            y=2.0 | nbody_system.length,
            z=3.0 | nbody_system.length,
            vx=1.0  | nbody_system.speed,
            vy=2.0 | nbody_system.speed,
            vz=3.0 | nbody_system.speed,
            mass=1.0 | nbody_system.mass,
            radius=4.0 | nbody_system.length,
        )
        instance.particles.add_particle(p) 
        self.assertEqual(instance.particles[0].radius, 0.0 | nbody_system.length)
        self.assertEqual(instance.particles[1].radius, 0.0 | nbody_system.length)
        self.assertEqual(instance.particles[2].radius, 4.0 | nbody_system.length)

        instance.stop()


    def test16(self):
        particles = new_plummer_model(200)
        particles.scale_to_standard()
        instance = Hermite()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000 | nbody_system.length**2
        instance.particles.add_particles(particles)

        x = numpy.arange(-1,1,0.1) | nbody_system.length
        zero = numpy.zeros(len(x)) | nbody_system.length
        potential0 =  instance.get_potential_at_point(zero, x , zero, zero)
        instance.stop()
        for n in (2, 3, 4):

            instance = Hermite(number_of_workers=n)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 0.00000 | nbody_system.length**2
            instance.particles.add_particles(particles)
            potential =  instance.get_potential_at_point(zero, x , zero, zero)

            self.assertAlmostRelativeEquals(potential0, potential, 8)
            instance.stop()



    def test17(self):
        particles = new_plummer_model(50)
        particles.scale_to_standard()
        instance = Hermite()
        instance.parameters.epsilon_squared = 0.22000 | nbody_system.length**2
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 50)

        instance.reset()

        self.assertAlmostRelativeEquals(instance.parameters.epsilon_squared , 0.22000 | nbody_system.length**2)
        self.assertEqual(len(instance.particles), 0)
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 50)
        instance.stop()


    def test18(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,1.0] | nbody_system.length
        particles.y = 0.0 | nbody_system.length
        particles.z = 0.0 | nbody_system.length
        particles.vx =  0.0 | nbody_system.speed
        particles.vy =  0.0 | nbody_system.speed
        particles.vz =  0.0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        instance = Hermite()
        instance.particles.add_particles(particles) 
        instance.commit_particles()
        self.assertEqual(instance.particles[0].radius, 0.0 | nbody_system.length)
        instance.parameters.end_time_accuracy_factor = 1.0
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostRelativeEquals(instance.model_time, 0.10563767746 |nbody_system.time, 5)
        instance.parameters.end_time_accuracy_factor = -1.0
        instance.evolve_model(0.3 | nbody_system.time)
        self.assertAlmostRelativeEquals(instance.model_time, 0.266758127609 |nbody_system.time, 5)
        instance.parameters.end_time_accuracy_factor = 0.0
        instance.evolve_model(0.4 | nbody_system.time)
        self.assertAlmostRelativeEquals(instance.model_time, 0.4 |nbody_system.time, 6)
        instance.parameters.end_time_accuracy_factor = -0.5
        instance.evolve_model(0.5 | nbody_system.time)
        self.assertAlmostRelativeEquals(instance.model_time, 0.48974930698 |nbody_system.time, 6)
        instance.parameters.end_time_accuracy_factor = +0.5
        instance.evolve_model(0.6 | nbody_system.time)
        self.assertAlmostRelativeEquals(instance.model_time, 0.6042733579 |nbody_system.time, 6)

        instance.stop()


    def test19(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,200.0] | nbody_system.length
        particles.y = 0.0 | nbody_system.length
        particles.z = 0.0 | nbody_system.length
        particles.vx =  0.0 | nbody_system.speed
        particles.vy =  0.0 | nbody_system.speed
        particles.vz =  0.0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        instance = Hermite()
        instance.particles.add_particles(particles) 
        instance.commit_particles()
        self.assertEqual(instance.particles[0].radius, 0.0 | nbody_system.length)
        instance.parameters.end_time_accuracy_factor = 0.0
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostRelativeEquals(instance.model_time, 0.1 |nbody_system.time, 5)

        instance.stop()

    def test20(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = Hermite(convert_nbody)
        hermite.initialize_code()
        hermite.parameters.epsilon_squared = 0.0 | units.AU**2
        hermite.parameters.end_time_accuracy_factor = 0.0
        hermite.parameters.is_time_reversed_allowed = True

        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        hermite.particles.add_particles(stars)

        hermite.evolve_model(365.0 | units.day)
        hermite.particles.copy_values_of_all_attributes_to(stars)

        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostRelativeEqual(position_at_start, position_after_full_rotation, 6)

        hermite.evolve_model(365.0 - (365.0 / 2) | units.day)

        hermite.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation_backward = earth.position.value_in(units.AU)[0]
        self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation_backward, 4)

        hermite.evolve_model(365.0 | units.day)

        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostRelativeEqual(position_at_start, position_after_full_rotation, 6)

        hermite.cleanup_code()

        hermite.stop()

    def test21(self):
        import pickle
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = Hermite(convert_nbody)
        encoded_interface = pickle.dumps(hermite,0)
        decoded_interface = pickle.loads(encoded_interface)


    def test22(self):
        hermite = Hermite()
        hermite.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = datamodel.Particles(2)
        particles.position = ([0,0,0], [1,0,0] )| nbody_system.length
        particles.velocity = ([-2,0,0], [2,0,0] )| nbody_system.speed
        particles.radius = 0| nbody_system.length
        particles.mass = 0.1| nbody_system.mass

        hermite.particles.add_particles(particles)
        hermite.stopping_conditions.out_of_box_detection.enable()
        hermite.parameters.stopping_conditions_out_of_box_size = 2 | nbody_system.length
        hermite.parameters.stopping_conditions_out_of_box_use_center_of_mass = False
        hermite.evolve_model(1 | nbody_system.time)
        print(hermite.particles.x)
        print(hermite.particles.key, particles[1].key)
        print(hermite.stopping_conditions.out_of_box_detection.particles(0))
        self.assertTrue(hermite.stopping_conditions.out_of_box_detection.is_set())
        self.assertEqual(len(hermite.stopping_conditions.out_of_box_detection.particles(0)), 1)
        self.assertEqual(hermite.stopping_conditions.out_of_box_detection.particles(0)[0].key, particles[1].key)
        hermite.stop()


    def test23(self):
        hermite = Hermite()
        hermite.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = datamodel.Particles(1)
        particles.position = ([0,0,0] )| nbody_system.length
        particles.velocity = ([1,0,0] )| nbody_system.speed
        particles.radius = 0| nbody_system.length
        particles.mass = 0.1| nbody_system.mass

        hermite.particles.add_particles(particles)
        hermite.evolve_model(1 | nbody_system.time)
        print(hermite.particles.x)
        self.assertAlmostRelativeEquals(hermite.model_time, 1 | nbody_system.time)
        self.assertAlmostRelativeEquals(hermite.particles[0].x, 1 | nbody_system.length)
        hermite.evolve_model(1.5 | nbody_system.time)
        print(hermite.particles.x)
        self.assertAlmostRelativeEquals(hermite.model_time, 1.5 | nbody_system.time)
        self.assertAlmostRelativeEquals(hermite.particles[0].x, 1.5 | nbody_system.length)
        hermite.stop()


    def test24(self):
        hermite = Hermite(reuse_worker=True)
        channel1 = hermite.legacy_interface.channel
        hermite.stop()
        hermite = Hermite(reuse_worker=True)
        channel2 = hermite.legacy_interface.channel
        hermite.stop()
        self.assertEqual(id(channel1), id(channel2))

    def test25(self):
        hermite = Hermite()
        hermite.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        particles = datamodel.Particles(10)
        particles.position = ([0,0,0] )| nbody_system.length
        particles.velocity = ([1,0,0] )| nbody_system.speed
        particles.radius = 0| nbody_system.length
        particles.mass = 0.1| nbody_system.mass
        particles.x = numpy.linspace(1, 10, 10) | nbody_system.length 
        particles.vx = numpy.linspace(1, 5, 10) | nbody_system.speed 

        hermite.particles.add_particles(particles)

        request = hermite.particles.get_values_in_store_async(None, ["x"])
        request.wait()
        print(request.result())
        self.assertEqual(request.result()[0], particles.x)
        request = hermite.particles.get_values_in_store_async(None, ["x", "vx"])
        request.wait()
        print(request.result())
        self.assertEqual(request.result()[0], particles.x)
        self.assertEqual(request.result()[1], particles.vx)
        p = particles.copy()
        channel = hermite.particles.new_channel_to(p)
        p.x = 0 | nbody_system.length
        p.vx = 0 | nbody_system.speed
        request = channel.copy_attributes_async(("x","vx",), async_get=True)
        request.wait()
        self.assertEqual(p.x, particles.x)
        self.assertEqual(p.vx, particles.vx)
        p.x = 0 | nbody_system.length
        p.vx = 0 | nbody_system.speed
        channel = p.new_channel_to(hermite.particles)
        request = channel.copy_attributes_async(("x", "y", "z","vx","vy","vz"), async_get=False, async_set=True)
        request.wait()
        self.assertEqual(p.x, hermite.particles.x)
        self.assertEqual(p.vx, hermite.particles.vx)
        channel = p.new_channel_to(particles)
        request = channel.copy_attributes_async(("x", "y", "z","vx","vy","vz"), async_get=False, async_set=True)
        request.wait()
        self.assertEqual(p.x, particles.x)
        self.assertEqual(p.vx, particles.vx)
        request = channel.copy_attributes_async(("x", "y", "z","vx","vy","vz"), async_get=True, async_set=False)
        request.wait()
        self.assertEqual(p.x, particles.x)
        self.assertEqual(p.vx, particles.vx)














