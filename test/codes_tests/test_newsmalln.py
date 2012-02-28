from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.newsmallN.interface import SmallNInterface, SmallN

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.ic import plummer
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestSmallNInterface(TestWithMPI):

    def test0(self):
        instance = SmallNInterface()
        instance.stop()
    
    def test1(self):
        instance = SmallNInterface()
        instance.initialize_code()
    
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEquals(1, res1['index_of_the_particle'])
        self.assertEquals(2, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
    
        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals(0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])
        self.assertEquals(2.0,  retrieved_state1['radius'])
        self.assertEquals(5.0,  retrieved_state2['radius'])
    
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        instance = SmallNInterface()
        instance.initialize_code()
        self.skip("index of the next particle not implemented correctly yet")
        for i in [0, 1, 2]:
            temp_particle = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEquals(i+1, temp_particle['index_of_the_particle'])
            
        instance.delete_particle(2)
      
        self.assertEquals(2, instance.get_number_of_particles()['number_of_particles'])
        
        self.assertEquals(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEquals(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle'])
        self.assertEquals(0, instance.get_index_of_next_particle(1)['__result'])
        self.assertEquals(-1, instance.get_index_of_next_particle(3)['__result'])
        self.assertEquals(1, instance.get_index_of_next_particle(2)['__result'])
        instance.stop()
        

    def test5(self):
        smalln = SmallNInterface()
        smalln.initialize_code()
        
        smalln.new_particle([10,20],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0],[1,1])
        retrieved_state = smalln.get_state(1)
        
        self.assertEquals(10.0,  retrieved_state['mass'])
        self.assertEquals(1, retrieved_state['radius'])
    
        retrieved_state = smalln.get_state([1,2])
        self.assertEquals(20.0,  retrieved_state['mass'][1])
        self.assertEquals(smalln.get_number_of_particles()['number_of_particles'], 2)
        smalln.cleanup_code() 
        smalln.stop()

    def test6(self):
        smalln = SmallNInterface()
        smalln.initialize_code()
        
        smalln.new_particle([10,10],[-1,1],[0,0], [0,0], [0,0], [0,0], [0,0], [1,1])
        retrieved_state = smalln.get_state(1)
        
        retr = smalln.get_potential_at_point(0.01, 0, 0, 0)
        self.assertEqual(retr['__result'], -1)
        smalln.cleanup_code()
        smalln.stop()
       
    def xtest7(self):
        instance = SmallNInterface()
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
        instance.stop()
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 10.0]) / 2.0)
        

    def xtest8(self):
        instance = SmallNInterface()
        instance.initialize_code()
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 1.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -1.0 / numpy.sqrt(2.0**2), 8)
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        instance.stop()
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 1.0]) / 2.0)
        
    
    def test9(self):
        print "Test SmallNInterface evolve_model"
        instance = SmallNInterface()
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_eta(0.001))
        self.assertEquals(0, instance.commit_parameters())
        
        # Set up an equal-mass binary on a circular orbit:
        self.assertEquals([1, 0], instance.new_particle(0.5,  0.5, 0, 0,  0, 0.5, 0, 0.01).values())
        self.assertEquals([2, 0], instance.new_particle(0.5,  -0.5, 0, 0,  0,-0.5, 0, 0.01).values())
        self.assertEquals(0, instance.commit_particles())
        
        self.assertEquals(0, instance.evolve_model(math.pi))
        for result, expected in zip(instance.get_position(1).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 3)
        for result, expected in zip(instance.get_position(2).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 3)
        
        self.assertEquals(0, instance.evolve_model(2 * math.pi))
        #print instance.get_time()
        #print instance.get_position(1), instance.get_velocity(1)
        #print instance.get_position(2)
        #for result, expected in zip(instance.get_position(1).values(), [0.5, 0.0, 0.0, 0]):
        #    self.assertAlmostEquals(result, expected, 3)
        #for result, expected in zip(instance.get_position(2).values(), [-0.5, 0.0, 0.0, 0]):
        #    self.assertAlmostEquals(result, expected, 3)
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    


class TestSmallN(TestWithMPI):
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
        self.skip("newsmalln starts every evolve at time 0")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        smalln = SmallN(convert_nbody)
        smalln.initialize_code()
        smalln.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        smalln.particles.add_particles(stars)
        
        smalln.evolve_model(365.0 | units.day)
        smalln.particles.copy_values_of_all_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        smalln.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        smalln.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        smalln.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        smalln.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        smalln.cleanup_code()
        
        smalln.stop()
        

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        instance = SmallN(convert_nbody)
        instance.initialize_code()
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.particles.add_particles(stars)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_all_attributes_to(stars)
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
            output_file = os.path.join(test_results_path, "smalln-earth-sun2.svg")
            figure.savefig(output_file)
        
        
        
        instance.cleanup_code()
        instance.stop()

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = SmallN(convert_nbody)
        instance.initialize_code()
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
        
        instance.stop()
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = SmallN(convert_nbody)
        instance.initialize_code()
        
        particles = datamodel.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        
        self.assertEquals(instance.get_mass(1), 17.0| units.kg) 
        self.assertEquals(instance.get_mass(2), 33.0| units.kg)  
        instance.stop()

    def test5(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = SmallN(convert_nbody)
        instance.initialize_code()
        
        particles = datamodel.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
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
        
        self.assertEquals(curr_state[0], 16|units.kg, 8)
    
    def test6(self):
        print "Test6: Testing SmallN parameters"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = SmallN(convert_nbody)
        instance.initialize_code()
        
       
        value = instance.get_eta()
        self.assertEquals(0.14 | units.none, value)
        self.assertAlmostEquals(0.14 | units.none, instance.parameters.timestep_parameter, in_units=units.none)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.timestep_parameter = x
            self.assertAlmostEquals(x | units.none, instance.parameters.timestep_parameter, in_units=units.none)
        
        
        value = instance.get_time()
        self.assertEquals(0| units.yr, value)
        
        value = instance.get_gamma()
        self.assertEquals(1e-6| units.none, value)
        self.assertAlmostEquals(1e-6| units.none, instance.parameters.unperturbed_threshold, in_units=units.none)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.unperturbed_threshold = x
            self.assertAlmostEquals(x | units.none, instance.parameters.unperturbed_threshold, in_units=units.none)
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

        instance = SmallN()
        instance.particles.add_particles(particles) 
        instance.commit_particles()
        self.assertEquals(instance.particles[0].radius, 0.0 | nbody_system.length)
        p = datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 3.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
            radius = 4.0 | nbody_system.length,
        )
        instance.particles.add_particle(p) 
        self.assertEquals(instance.particles[0].radius, 0.0 | nbody_system.length)
        self.assertEquals(instance.particles[1].radius, 0.0 | nbody_system.length)
        self.assertEquals(instance.particles[2].radius, 4.0 | nbody_system.length)
        
        instance.stop()
