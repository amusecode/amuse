from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.hermite0.interface import HermiteInterface, Hermite

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.ext import plummer

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestHermiteInterface(TestWithMPI):

    def test0(self):
        instance = HermiteInterface()
        self.assertTrue("Hut" in instance.all_literature_references_string())
    
    def test1(self):
        instance = HermiteInterface()
        instance.initialize_code()
    
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEquals(0, res1['index_of_the_particle'])
        self.assertEquals(1, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)
    
        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals(0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])
    
        instance.cleanup_code()
        del instance

    def test2(self):
        instance = HermiteInterface()
        instance.initialize_code()

        for i in [0, 1, 2]:
            temp_particle = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEquals(i, temp_particle['index_of_the_particle'])
            
        instance.delete_particle(1)
      
        self.assertEquals(2, instance.get_number_of_particles()['number_of_particles'])
        
        self.assertEquals(0, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEquals(2, instance.get_index_of_next_particle(0)['index_of_the_next_particle'])
        self.assertEquals(0, instance.get_index_of_next_particle(0)['__result'])
        self.assertEquals(-1, instance.get_index_of_next_particle(1)['__result'])
        self.assertEquals(1, instance.get_index_of_next_particle(2)['__result'])
        
    def test3(self):
        hermite = HermiteInterface()
        hermite.eps2 = 0.101
        self.assertEquals(0.101, hermite.eps2)
        hermite.eps2 = 0.110
        self.assertEquals(0.110, hermite.eps2)
        hermite.stop()

    def test4(self):
        hermite = HermiteInterface()
        hermite.flag_collision = 1
        self.assertEquals(1, hermite.flag_collision)
        hermite.flag_collision = 0
        self.assertEquals(0, hermite.flag_collision)
        hermite.stop()

    def test5(self):
        hermite = HermiteInterface()
        hermite.initialize_code()
        
        hermite.new_particle([10,20],[1,1],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0])
        retrieved_state = hermite.get_state(0)
        
        self.assertEquals(10.0,  retrieved_state['mass'])
        self.assertEquals(1, retrieved_state['radius'])
    
        retrieved_state = hermite.get_state([0,1])
        self.assertEquals(20.0,  retrieved_state['mass'][1])
        self.assertEquals(hermite.get_number_of_particles()['number_of_particles'], 2)
        hermite.cleanup_code() 

    def test6(self):
        hermite = HermiteInterface()
        hermite.initialize_code()
        
        hermite.new_particle([10,10],[1,1],[-1,1],[0,0], [0,0], [0,0], [0,0], [0,0])
        retrieved_state = hermite.get_state(0)
        
        retr = hermite.get_potential_at_point(0.01, 0,0,0)
        self.assertEqual(retr['phi'], -20.0)
        hermite.cleanup_code()
       
    def test7(self):
        instance = HermiteInterface()
        instance.initialize_code()
        instance.set_eps2(0.1 * 0.1)
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2 + 0.1**2), 8)
        
        potential, errorcode = instance.get_potential(id2)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2 + 0.1**2), 8)
        
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 10.0]) / 2.0)
        

    def test8(self):
        instance = HermiteInterface()
        instance.initialize_code()
        instance.set_eps2(0)
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 1.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -1.0 / numpy.sqrt(2.0**2), 8)
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 1.0]) / 2.0)
        
class TestHermite(TestWithMPI):
    def new_system_of_sun_and_earth(self):
        stars = core.Stars(2)
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
        hermite.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        hermite.particles.add_particles(stars)
        
        hermite.evolve_model(365.0 | units.day)
        hermite.particles.copy_values_of_state_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        hermite.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        hermite.particles.copy_values_of_state_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        hermite.particles.copy_values_of_state_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
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
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_state_attributes_to(stars)
            stars.savepoint()
        
        if HAS_MATPLOTLIB:
            figure = pyplot.figure()
            plot = figure.add_subplot(1,1,1)
            
            x_points = earth.get_timeline_of_attribute("x")
            y_points = earth.get_timeline_of_attribute("y")
            
            x_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), x_points)
            y_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), y_points)
            
            plot.scatter(x_points_in_AU,y_points_in_AU, color = "b", marker = 'o')
            
            plot.set_xlim(-1.5, 1.5)
            plot.set_ylim(-1.5, 1.5)
               
            
            test_results_path = self.get_path_to_results()
            output_file = os.path.join(test_results_path, "hermite-earth-sun2.svg")
            figure.savefig(output_file)
        
        
        
        instance.cleanup_code()
        del instance

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Hermite(convert_nbody, redirection="none")
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | units.AU**2
        instance.dt_dia = 5000
        
        stars = core.Stars(2)
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
            instance.particles.copy_values_of_state_attributes_to(stars)
            stars.savepoint()
        
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Hermite(convert_nbody, redirection="none")
        instance.initialize_code()
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        
        self.assertEquals(instance.get_mass(0), 17.0| units.kg) 
        self.assertEquals(instance.get_mass(1), 33.0| units.kg)  

    def test5(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Hermite(convert_nbody)
        instance.initialize_code()
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        instance.set_state(1, 16|units.kg, 20.0|units.m, 
                           20.0|units.m, 40.0|units.m, 60.0|units.m, 
                           1.0|units.ms, 1.0|units.ms, 1.0|units.ms)
        
        curr_state =  instance.get_state(1)
        
        self.assertEquals(curr_state[0], 16|units.kg, 8)
    
    def test6(self):
        print "Test6: Testing Hermite parameters"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = Hermite(convert_nbody)
        
        (value, error) = instance.get_eps2()
        self.assertEquals(0, error)
        self.assertEquals(0.0, value)
        self.assertAlmostEquals(0.0 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEquals(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        
        (value, error) = instance.get_dt_param()
        self.assertEquals(0, error)
        self.assertEquals(0.03, value)
        self.assertAlmostEquals(0.03 | units.none, instance.parameters.dt_param, in_units=units.none)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.dt_param = x | units.none
            self.assertAlmostEquals(x | units.none, instance.parameters.dt_param, in_units=units.none)
        
        (value, error) = instance.get_dt_dia()
        self.assertEquals(0, error)
        self.assertEquals(1.0, value)
        self.assertAlmostEquals(1.0 | units.yr, instance.parameters.dt_dia, in_units=units.yr)
        for x in [0.1, 10.0, 100.0]:
            instance.parameters.dt_dia = x | units.yr
            self.assertAlmostEquals(x | units.yr, instance.parameters.dt_dia, in_units=units.yr)
        
        (value, error) = instance.get_time()
        self.assertEquals(0, error)
        self.assertEquals(0.0, value)
        self.assertAlmostEquals(0.0 | units.yr, instance.parameters.time, in_units=units.yr)
        for x in [1.0, 10.0, 100.0]:
            instance.parameters.time = x | units.yr
            self.assertAlmostEquals(x | units.yr, instance.parameters.time, in_units=units.yr)
        instance.stop()
    
    def test7(self):
        print "Test7: Testing effect of Hermite parameter epsilon_squared"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        
        particles = core.Particles(2)
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
        
        initial_direction = math.atan((earth.velocity[0]/earth.velocity[1]).value_in(units.none))
        final_direction = []
        for log_eps2 in range(-9,10,2):
            instance = Hermite(convert_nbody)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            final_direction.append(math.atan((instance.particles[1].velocity[0]/
                instance.particles[1].velocity[1]).value_in(units.none)))
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEquals(abs(final_direction[0]), abs(initial_direction+math.pi/2.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEquals(final_direction[-1], initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(final_direction[i+1]-final_direction[i]) for i in range(len(final_direction)-1)]
        self.assertEquals(delta[len(final_direction)/2 -1], max(delta))
        
    def test8(self):
        particles = core.Particles(2)
        particles.x = [
            0.0,1.0, 
            #5,7,
            #10,12,
            #15,17,
            #20,22
        ] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.75 | nbody_system.length
        particles.vx =  0.1 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 0 | nbody_system.mass
       
        instance = Hermite()
        instance.initialize_code()
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.collision_detection.enable()
        instance.evolve_model(0.1 | nbody_system.time)
        print instance.model_time
        self.assertTrue(instance.stopping_conditions.collision_detection.is_set())
        print instance.stopping_conditions.collision_detection.particles(0).key
        print instance.stopping_conditions.collision_detection.particles(1).key
        self.assertEquals(len(instance.stopping_conditions.collision_detection.particles(0)), 1 )
        p0 =  instance.stopping_conditions.collision_detection.particles(0)[0]
        p1 =  instance.stopping_conditions.collision_detection.particles(1)[0]
        self.assertNotEquals(p0, p1)
        print p0.x, p1.x
        self.assertTrue(p1.x - p0.x < 1.5| nbody_system.length)
        instance.stop()
        
    def test9(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        instance = Hermite(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.pair_factor = 1e14 | units.none
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        instance.particles.add_particles(stars)
        instance.stopping_conditions.pair_detection.enable()
        instance.evolve_model(365.0 | units.day)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        instance.particles.copy_values_of_state_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
    
        print instance.model_time
        
        print instance.stopping_conditions.pair_detection.particles(0).key
        print instance.stopping_conditions.pair_detection.particles(1).key
        
        instance.cleanup_code()
        
        instance.stop()
        
    def test10(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        instance = Hermite(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.stopping_conditions_number_of_steps = 10
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps,10|units.none)
    
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        instance.particles.add_particles(stars)
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(365.0 | units.day)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        instance.particles.copy_values_of_state_attributes_to(stars)
        
        instance.cleanup_code()
        
        instance.stop()

    def test11(self):
        particles = core.Particles(2)
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
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 2|units.none)
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody_system.time)
        
        instance.stop()

    def test12(self):
        particles = core.Particles(2)
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
        self.assertEquals(instance.parameters.stopping_conditions_out_of_box_size, .5 | nbody_system.length)
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.out_of_box_detection.enable()
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.out_of_box_detection.is_set())
        self.assertAlmostEqual(instance.stopping_conditions.out_of_box_detection.particles(0).x, 1.0 |nbody_system.length, 3)
        instance.stop()

    def test13(self):
        particles = plummer.new_plummer_sphere(31)
       
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
            self.assertAlmostRelativeEqual(expected_positions, positions_per_workers[index])
    
    
