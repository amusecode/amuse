from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.smalln.interface import SmallNInterface, SmallN

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse import io
from amuse.datamodel import trees
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
        
        self.assertEqual(1, res1['index_of_the_particle'])
        self.assertEqual(2, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
    
        self.assertEqual(11.0,  retrieved_state1['mass'])
        self.assertEqual(21.0,  retrieved_state2['mass'])
        self.assertEqual(0.0,  retrieved_state1['x'])
        self.assertEqual(10.0,  retrieved_state2['x'])
        self.assertEqual(2.0,  retrieved_state1['radius'])
        self.assertEqual(5.0,  retrieved_state2['radius'])
    
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        instance = SmallNInterface()
        instance.initialize_code()
        self.skip("index of the next particle not implemented correctly yet")
        for i in [0, 1, 2]:
            temp_particle = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEqual(i+1, temp_particle['index_of_the_particle'])
            
        instance.delete_particle(2)
      
        self.assertEqual(2, instance.get_number_of_particles()['number_of_particles'])
        
        self.assertEqual(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEqual(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle'])
        self.assertEqual(0, instance.get_index_of_next_particle(1)['__result'])
        self.assertEqual(-1, instance.get_index_of_next_particle(3)['__result'])
        self.assertEqual(1, instance.get_index_of_next_particle(2)['__result'])
        instance.cleanup_code()
        instance.stop()
        

    def test5(self):
        smalln = SmallNInterface()
        smalln.initialize_code()
        
        smalln.new_particle([10,20],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0],[1,1])
        retrieved_state = smalln.get_state(1)
        
        self.assertEqual(10.0,  retrieved_state['mass'])
        self.assertEqual(1, retrieved_state['radius'])
    
        retrieved_state = smalln.get_state([1,2])
        self.assertEqual(20.0,  retrieved_state['mass'][1])
        self.assertEqual(smalln.get_number_of_particles()['number_of_particles'], 2)
        smalln.cleanup_code() 
        smalln.stop()

    def test6(self):
        smalln = SmallNInterface()
        smalln.initialize_code()
        
        smalln.new_particle([10,10],[-1,1],[0,0], [0,0], [0,0], [0,0], [0,0], [1,1])
        retrieved_state = smalln.get_state(1)
        
        self.assertFalse(hasattr(smalln, 'get_potential_at_point'))
        #retr = smalln.get_potential_at_point(0.01, 0, 0, 0)
        #self.assertEqual(retr['__result'], -1)
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
        

    def xtest8(self):
        instance = SmallNInterface()
        instance.initialize_code()
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 1.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
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
        print("Test SmallNInterface evolve_model")
        instance = SmallNInterface()
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.set_eta(0.001))
        self.assertEqual(0, instance.commit_parameters())
        
        # Set up an equal-mass binary on a circular orbit:
        self.assertEqual([1, 0], list(instance.new_particle(0.5,  0.5, 0, 0,  0, 0.5, 0, 0.01).values()))
        self.assertEqual([2, 0], list(instance.new_particle(0.5,  -0.5, 0, 0,  0,-0.5, 0, 0.01).values()))
        self.assertEqual(0, instance.commit_particles())
        
        self.assertEqual(0, instance.evolve_model(math.pi))
        for result, expected in zip(instance.get_position(1).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(instance.get_position(2).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        
        self.assertEqual(0, instance.evolve_model(2 * math.pi))
        #print instance.get_time()
        #print instance.get_position(1), instance.get_velocity(1)
        #print instance.get_position(2)
        #for result, expected in zip(instance.get_position(1).values(), [0.5, 0.0, 0.0, 0]):
        #    self.assertAlmostEquals(result, expected, 3)
        #for result, expected in zip(instance.get_position(2).values(), [-0.5, 0.0, 0.0, 0]):
        #    self.assertAlmostEquals(result, expected, 3)
        
        self.assertEqual(0, instance.cleanup_code())
        instance.cleanup_code()
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
            
            x_points_in_AU = [t_x[1].value_in(units.AU) for t_x in x_points]
            y_points_in_AU = [t_x1[1].value_in(units.AU) for t_x1 in y_points]
            
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
        self.assertEqual(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 2)
        
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        
        self.assertEqual(instance.get_mass(1), 17.0| units.kg) 
        self.assertEqual(instance.get_mass(2), 33.0| units.kg)  
        instance.stop()

    def test5(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = SmallN(convert_nbody)
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
        print("Test6: Testing SmallN parameters")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = SmallN(convert_nbody)
        instance.initialize_code()
        
       
        value = instance.get_eta()
        self.assertEqual(0.14, value)
        self.assertAlmostEqual(0.14, instance.parameters.timestep_parameter)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.timestep_parameter = x
            self.assertAlmostEqual(x, instance.parameters.timestep_parameter)
        
        
        value = instance.get_time()
        self.assertEqual(0| units.yr, value)
        
        value = instance.get_gamma()
        self.assertEqual(1e-6, value)
        self.assertAlmostEqual(1e-6, instance.parameters.unperturbed_threshold)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.unperturbed_threshold = x
            self.assertAlmostEqual(x, instance.parameters.unperturbed_threshold)
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
        self.assertEqual(instance.particles[0].radius, 0.0 | nbody_system.length)
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
        self.assertEqual(instance.particles[0].radius, 0.0 | nbody_system.length)
        self.assertEqual(instance.particles[1].radius, 0.0 | nbody_system.length)
        self.assertEqual(instance.particles[2].radius, 4.0 | nbody_system.length)
        
        instance.stop()

    def assertEarthAndMoonWasDetectedAsBinary(self, set, stars):
        alltrees = trees.BinaryTreesOnAParticleSet(set, 'child1', 'child2')
        roots = list(alltrees.iter_roots())
        self.assertEqual(len(roots), 1)
        root = roots[0].particle
        earth_and_moon = None
        if not root.child1.child1 is None:
            earth_and_moon = root.child1
        else:
            earth_and_moon = root.child2
        
        self.assertAlmostRelativeEquals(earth_and_moon.mass, stars[1].mass + stars[2].mass)
        
        
    def test16(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        smalln = SmallN(convert_nbody)
        smalln.initialize_code()
        smalln.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        
        moon = datamodel.Particle()
        moon.mass = units.kg(7.3477e22)
        moon.radius = units.km(1737.10) 
        moon.position = units.km(numpy.array((149.5e6 + 384.399 ,0.0,0.0)))
        moon.velocity = units.ms(numpy.array((0.0,29800 + 1022,0.0)))
        
        stars.add_particle(moon)
        
        earth = stars[1]
                
        smalln.particles.add_particles(stars)
        
        smalln.evolve_model(365.0 | units.day)
        smalln.update_particle_tree()
        smalln.update_particle_set()
        

        self.assertEqual(len(smalln.particles), 5)
        
        self.assertEarthAndMoonWasDetectedAsBinary(smalln.particles, stars)
        
        
        inmemory = smalln.particles.copy()
        self.assertEarthAndMoonWasDetectedAsBinary(inmemory, stars)
        
        test_results_path = self.get_path_to_results()
        output_file = os.path.join(test_results_path, "newsmalln-test16.hdf5")
        if os.path.exists(output_file):
            os.remove(output_file)
            
        io.write_set_to_file(smalln.particles, output_file, "hdf5")
        fromfile = io.read_set_from_file(output_file, "hdf5")
        self.assertEarthAndMoonWasDetectedAsBinary(fromfile, stars)
        smalln.stop()

        
   
    def test17(self):
        
        particles = datamodel.Particles(keys=[1,2,3,4,5,6,7])
        particles.mass = 0.001 | nbody_system.mass
        particles.radius = 0.1 | nbody_system.length
        particles.x = [
            -100.5, -99.5, 
              -0.5, 0.5, 
              99.5, 100.5,
             120.0
        ] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [
            [2, 0, 0], [-2, 0, 0],
            [2, 0, 0], [-2, 0, 0],
            [2, 0, 0], [-2, 0, 0],
            [-4, 0, 0]
        ] | nbody_system.speed
        
        instance = SmallN()
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
        instance.stop()
        
    
    def test18(self):
        
        particles = datamodel.Particles(keys=[1,2])
        particles.mass = 1 | nbody_system.mass
        particles.radius = 0.1 | nbody_system.length
        particles.x = [1, -1] | nbody_system.length
        particles.y = [1, -1] | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[-1, 0, 0], [1, 0, 0]] | nbody_system.speed
       
        instance = SmallN()
        instance.particles.add_particles(particles)
        
        stopping_condition = instance.stopping_conditions.interaction_over_detection
        stopping_condition.enable()
        
        
        instance.evolve_model(10.0 | nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        self.assertTrue(instance.model_time < 11.0 | nbody_system.time)
        instance.stop() 
        

    def test19(self):
        
        particles = datamodel.Particles(keys=[1,2, 3])
        particles.mass = 1 | nbody_system.mass
        particles.radius = 0.1 | nbody_system.length
        particles.x = [1, -1, 0] | nbody_system.length
        particles.y = [1, -1, 0] | nbody_system.length
        particles.z = [0, 0, 1]| nbody_system.length
        particles.velocity = [[-1, 0, 0], [1, 0, 0],[0,0,-10]] | nbody_system.speed
       
        instance = SmallN()
        instance.particles.add_particles(particles)
        
        stopping_condition = instance.stopping_conditions.interaction_over_detection
        stopping_condition.enable()
        
        
        instance.evolve_model(10.0 | nbody_system.time)
        self.assertTrue(stopping_condition.is_set())
        
        self.assertTrue(instance.model_time < 10.0 | nbody_system.time)
        instance.stop() 
    
    def test20(self):
        p = datamodel.Particles(3)

        p[0].mass   = 6.667e-01 | nbody_system.mass
        p[0].radius = 4.000e-03 | nbody_system.length
        p[0].x  = -1.309e+01 | nbody_system.length
        p[0].y  =  1.940e+01 | nbody_system.length
        p[0].z  = -1.163e+01 | nbody_system.length
        p[0].vx =  2.366e-01 | nbody_system.speed
        p[0].vy = -3.813e-01 | nbody_system.speed
        p[0].vz =  2.486e-01 | nbody_system.speed

        p[1].mass   = 3.333e-01 | nbody_system.mass
        p[1].radius = 1.000e-03 | nbody_system.length  
        p[1].x  = -1.506e+01 | nbody_system.length
        p[1].y  =  1.937e+01 | nbody_system.length
        p[1].z  = -1.163e+01 | nbody_system.length
        p[1].vx =  3.483e-01 | nbody_system.speed
        p[1].vy = -4.513e-01 | nbody_system.speed
        p[1].vz =  2.486e-01 | nbody_system.speed

        p[2].mass   = 5.000e-01 | nbody_system.mass
        p[2].radius = 2.000e-03 | nbody_system.length 
        p[2].x  =  2.749e+01 | nbody_system.length
        p[2].y  = -3.877e+01 | nbody_system.length
        p[2].z  =  2.325e+01 | nbody_system.length
        p[2].vx = -5.476e-01 | nbody_system.speed
        p[2].vy =  8.092e-01 | nbody_system.speed
        p[2].vz = -4.972e-01 | nbody_system.speed

        instance = SmallN()
        instance.initialize_code()
        instance.parameters.set_defaults
        N = 3
        t_begin = 0.0 | nbody_system.time
        t_end = 100.0 | nbody_system.time

        particles = p

        instance.particles.add_particles(particles)
        instance.commit_particles()


        sc = instance.stopping_conditions.collision_detection
        sc.enable()


        isCollision = False
        instance.evolve_model(t_end)
        isCollision = sc.is_set()
        
        instance.stop()
        self.assertTrue(isCollision, "no collision detected")
