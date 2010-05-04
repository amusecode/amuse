# -*- coding: utf-8 -*-
from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree
from amuse.support.data import core
from amuse.support.data import values
from amuse.support.data import particle_attributes
from amuse.support.units import constants
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.ext import plummer

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestMPIInterface(TestWithMPI):

    def test0(self):
        instance = BHTreeInterface()
        self.assertTrue("Barnes" in instance.all_literature_references_string())
        
    def test1(self):
        instance = BHTreeInterface()
        instance.initialize_code()
        instance.commit_parameters()
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        self.assertEquals(1, res1['index_of_the_particle'])
        self.assertEquals(2, res2['index_of_the_particle'])

        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)

        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals(0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])

        self.assertEquals(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        self.assertEquals(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle']) 
        
        instance.cleanup_code()
        del instance

    def test2(self):
        instance = BHTreeInterface()
        instance.initialize_code()

        instance.commit_parameters()
        for i in [1, 2, 3]:
            temp_particle = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEquals(i, temp_particle['index_of_the_particle'])
        
        instance.commit_particles()
        self.assertEquals(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        self.assertEquals(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle']) 
        self.assertEquals(3, instance.get_index_of_next_particle(2)['index_of_the_next_particle'])
            
        instance.delete_particle(1)
      
        self.assertEquals(2, instance.get_number_of_particles()['number_of_particles'])
        
        #the deletion does a swap, so 3 is copied to 1, (overwriting old 1 and treesize -> treesize-1
        self.assertEquals(3, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEquals(1, instance.get_index_of_next_particle(2)['__result'])

        instance.cleanup_code()
        del instance
        
    
        
    def test3(self):
        interface = BHTreeInterface()
        interface.eps2 = 0.101
        self.assertEquals(0.101, interface.eps2)
        interface.eps2 = 0.110
        self.assertEquals(0.110, interface.eps2)
        interface.stop()

    def test4(self):
        interface = BHTreeInterface()
        interface.flag_collision = 1
        self.assertEquals(1, interface.flag_collision)
        interface.flag_collision = 0
        self.assertEquals(0, interface.flag_collision)
        interface.stop()

    def test5(self):
        interface = BHTreeInterface()
        interface.initialize_code()
        
        interface.commit_parameters()
        interface.new_particle([10,20],[1,1],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0])
        interface.commit_particles()
        retrieved_state = interface.get_state(1)
        
        self.assertEquals(10.0,  retrieved_state['mass'])
        self.assertEquals(1, retrieved_state['radius'])

        retrieved_state = interface.get_state([1,2])
        self.assertEquals(20.0,  retrieved_state['mass'][1])
        self.assertEquals(interface.get_number_of_particles()['number_of_particles'], 2)
        interface.cleanup_code()
        del interface
    
    
    def test6(self):
        instance = BHTreeInterface()
        instance.initialize_code()
        instance.commit_parameters()
        
        ids = []
        for i in [1, 2, 3]:
            id, error = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            ids.append(id)
            interface.stop()
        print ids
        
        instance.commit_particles()
        
            
        instance.delete_particle(ids[0])
        id, error = instance.new_particle(mass = 4, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        self.failIfEqual(id, ids[-1])
        
        instance.cleanup_code()
        del instance


class TestAmuseInterface(TestWithMPI):
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = BHTree(convert_nbody)
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        
        stars = core.Stars(2)
        
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = [0.0,0.0,0.0] | units.m
        sun.velocity = [0.0,0.0,0.0] | units.ms
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = [149.5e6, 0.0, 0.0] | units.km
        earth.velocity = [0.0, 29800, 0.0] | units.ms

        #instance.particles.add_particles(stars)
        instance.particles.add_particles(stars)
        
        postion_at_start = earth.position.value_in(units.AU)[0]
        
        instance.evolve_model(365.0 | units.day)
        instance.update_particles(stars)
        
        postion_after_full_rotation = earth.position.value_in(units.AU)[0]
       
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 2)
        
        
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
         
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 1)
        instance.cleanup_code()
        del instance
        
    def test2(self):
        #not completed 
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = BHTree(convert_nbody)
        #instance.dt_dia = 1
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        #instance.timestep = 0.0001
        #instance.use_self_gravity = 0
        instance.commit_parameters()
        
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

        instance.particles.add_particles(stars)
        instance.commit_particles()
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.update_particles(stars)
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
            output_file = os.path.join(test_results_path, "bhtree-earth-sun.svg")
            figure.savefig(output_file)    
        
        instance.cleanup_code()
        del instance

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = BHTree(convert_nbody)
        #instance.dt_dia = 1
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        #instance.timestep = 0.0001
        #instance.use_self_gravity = 0
        instance.commit_parameters()
        
        
        stars = core.Stars(2)
        star1 = stars[0]
        star2 = stars[1]

        star1.mass = units.MSun(1.0)
        star1.position = units.AU(numpy.array((-.10,0.0,0.0)))
        star1.velocity = units.AUd(numpy.array((0.0,0.0,0.0)))
        star1.radius = units.RSun(1.0)

        star2.mass = units.MSun(1.0)
        star2.position = units.AU(numpy.array((.10,0.0,0.0)))
        star2.velocity = units.AUd(numpy.array((0.0,0.0,0.0)))
        star2.radius = units.RSun(100.0)
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
    
        for x in range(1,200,1):
            instance.evolve_model(x | units.day)
            instance.update_particles(stars)
            #instance.get_indices_of_colliding_particles()
            #print stars[0].position-stars[1].position
            stars.savepoint()
            
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = BHTree(convert_nbody)
        instance.commit_parameters()
        
        index = instance.new_particle(
            15.0 | units.kg,
            10.0 | units.m,
            10.0 | units.m, 20.0 | units.m, 30.0 | units.m,
            #1.0 | units.m/units.s, 1.0 | units.m/units.s, 3.0 | units.m/units.s
            0.0 | units.m/units.s, 0.0 | units.m/units.s, 0.0 | units.m/units.s
        )
        instance.commit_particles()
        self.assertEquals(instance.get_mass(index), 15.0| units.kg)
        
    def test5(self):

        instance = BHTree()
        instance.commit_parameters()
        
        index = instance.new_particle(
            15.0 | nbody_system.mass,
            10.0 | nbody_system.length,
            10.0 | nbody_system.length, 20.0 | nbody_system.length, 30.0 | nbody_system.length,
            1.0 | nbody_system.speed, 1.0 | nbody_system.speed, 3.0 | nbody_system.speed
        )
        instance.commit_particles()
        self.assertEquals(instance.get_mass(index), 15.0| nbody_system.mass)
        
    
    def test6(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)
        instance = BHTree(convert_nbody)
        instance.commit_parameters()
        
        indices = instance.new_particle(
            [15.0, 30.0] | units.kg,
            [10.0, 20.0] | units.m,
            [10.0, 20.0] | units.m, [20.0, 40.0] | units.m, [30.0, 50.0] | units.m,
            #1.0 | units.m/units.s, 1.0 | units.m/units.s, 3.0 | units.m/units.s
            [0.0, 0.01] | units.m/units.s, [0.0, 0.01] | units.m/units.s, [0.0, 0.01] | units.m/units.s
        )
        instance.commit_particles()
        
        self.assertEquals(instance.get_mass(indices[0]), 15.0| units.kg)
        self.assertEquals(instance.get_mass(indices)[0], 15.0| units.kg)
        
        try:
            instance.get_mass([4,5])
            self.fail("Should raise error, invalid index")
        except Exception, ex:
            self.assertEquals(str(ex), "Error when calling 'get_mass' of a 'BHTree', errorcode is -1")

        
    
    def test7(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = BHTree(convert_nbody)
        instance.commit_parameters()
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        self.assertEquals(instance.get_mass(1), 15.0| units.kg)  
        self.assertEquals(instance.get_position(1)[2], 30.0| units.m)   
        
        self.assertEquals(len(instance.particles), 2)
        
        
        self.assertEquals(instance.particles.mass[1] , 30.0 | units.kg) 
        self.assertEquals(instance.particles.position[1][2] , 60.0 | units.m)   
        
    def test8(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = BHTree(convert_nbody)
        instance.commit_parameters()
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        
        self.assertEquals(instance.get_mass(1), 17.0| units.kg) 
        
    def test9(self):
        instance = BHTree()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | nbody_system.length**2
        
        particles = core.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, 1.0 | nbody_system.length, zero, zero)
        self.assertAlmostEqual(fx, 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration, 3)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            potential1 = instance.get_potential_at_point(zero, x1, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            fx1, fy1, fz1 = instance.get_gravity_at_point(zero, x1, zero, zero)
            
            self.assertAlmostEqual(fy0, 0.0 | nbody_system.acceleration, 3)
            self.assertAlmostEqual(fz0, 0.0 | nbody_system.acceleration, 3)
            self.assertAlmostEqual(fy1, 0.0 | nbody_system.acceleration, 3)
            self.assertAlmostEqual(fz1, 0.0 | nbody_system.acceleration, 3)
            
            self.assertAlmostEqual(fx0, -1.0 * fx1, 5)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0, 2)
            self.assertAlmostEqual(potential0, potential1, 5)
            
    def test10(self):
        instance = BHTree()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | nbody_system.length**2
        instance.commit_parameters()
        
        
        particles = core.Particles(6)
        particles.mass = 1.0 | nbody_system.mass
        particles.radius =   0.00001 | nbody_system.length
        particles.position = [[-1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,1.0,0.0] ,[0.0,0.0,-1.0],[0.0,0.0,1.0]] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, zero, zero, zero)
        self.assertAlmostEqual(fx, 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration, 3)
        
        
        for position in (0.25, 0.5, 0.75):
            p0 = position | nbody_system.length
            p1 = -position | nbody_system.length
            for i in range(3):
                args0 = [zero] * 4
                args1 = [zero] * 4
                args0[1 + i] = p0
                args1[1 + i] = p1
                f0 = instance.get_gravity_at_point(*args0)
                f1 = instance.get_gravity_at_point(*args1)
                
                for j in range(3):
                    if j != i:
                        self.assertAlmostEqual(f0[j], 0.0 | nbody_system.acceleration, 3)
                        self.assertAlmostEqual(f1[j], 0.0 | nbody_system.acceleration, 3)
                    else:
                        self.assertAlmostEqual(f0[j], -1.0 * f1[j], 5)
        
        instance.stop()
        del instance
        
    def test11(self):
       
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = BHTree(convert_nbody)
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        
        copyof =  instance.particles.copy()
        
        self.assertAlmostEqual(30 | units.kg, copyof[1].mass, 6)
        
        copyof[1].mass = 35 | units.kg
        
        copyof.copy_values_of_state_attributes_to(instance.particles)
        
        self.assertAlmostEqual(35 | units.kg, instance.particles[1].mass, 6)

    def test12(self):
       
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = BHTree(convert_nbody)
        instance.commit_parameters()
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        copyof =  instance.particles.copy()
        
        instance.set_state(1, 16|units.kg, 20.0|units.m, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms)
        
        curr_state =  instance.get_state(1)
        
        self.assertEquals(curr_state[0], 16|units.kg, 8)

    def test13(self):
       
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.m)
        
        instance = BHTree(convert_nbody)
        instance.commit_parameters()
        
        particles = core.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [30.0, 30.0] | units.kg
        particles.radius =  [1.0, 1.0] | units.m
        particles.position = [[-10.0, 0.0, 0.0], [10.0, 0.0, 0.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s
        
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        copyof =  instance.particles.copy()
        
        com = instance.get_center_of_mass_position()
        self.assertAlmostEqual(com[0], values.new_quantity(0.0, units.m), constants.precision)
    
    def test14(self):
        print "Test14: Testing BHTree parameters (I)"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = BHTree(convert_nbody)
        
        (value, error) = instance.get_epsilon_squared()
        self.assertEquals(0, error)
        self.assertEquals(0.125, value)
        self.assertAlmostEquals(0.125 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEquals(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        
        instance.commit_parameters()
        (value, error) = instance.get_time_step()
        self.assertEquals(0, error)
        self.assertEquals(0.015625, value)
        self.assertAlmostEquals(0.015625 | units.yr, instance.parameters.timestep, in_units=units.yr)
        for x in [0.001, 0.01, 0.1]:
            instance.parameters.timestep = x | units.yr
            self.assertAlmostEquals(x | units.yr, instance.parameters.timestep, in_units=units.yr)
        
        (value, error) = instance.get_theta_for_tree()
        self.assertEquals(0, error)
        self.assertEquals(0.75, value)
        self.assertEquals(0.75 | units.none, instance.parameters.opening_angle)
        for x in [0.2, 0.5, 0.7]:
            instance.parameters.opening_angle = x | units.none
            self.assertEquals(x | units.none, instance.parameters.opening_angle)
        
        (value, error) = instance.get_use_self_gravity()
        self.assertEquals(0, error)
        self.assertEquals(1, value)
        self.assertEquals(1 | units.none, instance.parameters.use_self_gravity)
        for x in [0, 1]:
            instance.parameters.use_self_gravity = x | units.none
            self.assertEquals(x | units.none, instance.parameters.use_self_gravity)
        
        (value, error) = instance.get_ncrit_for_tree()
        self.assertEquals(0, error)
        self.assertEquals(1024, value)
        self.assertEquals(1024 | units.none, instance.parameters.ncrit_for_tree)
        for x in [512, 2048, 4096]:
            instance.parameters.ncrit_for_tree = x | units.none
            self.assertEquals(x | units.none, instance.parameters.ncrit_for_tree)
        
        (value, error) = instance.get_dt_dia()
        self.assertEquals(0, error)
        self.assertEquals(1.0, value)
        self.assertAlmostEquals(1.0 | units.yr, instance.parameters.dt_dia, in_units=units.yr)
        for x in [0.1, 10.0, 100.0]:
            instance.parameters.dt_dia = x | units.yr
            self.assertAlmostEquals(x | units.yr, instance.parameters.dt_dia, in_units=units.yr)
    
    def test15(self):
        print "Test15: Testing effect of BHTree parameter epsilon_squared"
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
            instance = BHTree(convert_nbody)
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
        
    
    def test16(self):
        numpy.random.seed(0)
        number_of_stars = 2
        stars = plummer.MakePlummerModel(number_of_stars).result
        stars.radius = 0.00001 | nbody_system.length
        stars.scale_to_standard()
        
        instance = BHTree()
        instance.initialize_code()
        instance.parameters.epsilon_squared = (1.0 / 20.0 / (number_of_stars**0.33333) | nbody_system.length)**2
        instance.parameters.timestep = 0.004 | nbody_system.time
        instance.parameters.timestep = 0.00001 | nbody_system.time
        instance.commit_parameters()
        print instance.parameters.timestep
        instance.particles.add_particles(stars)
        instance.commit_particles()
        energy_total_t0 = instance.potential_energy + instance.kinetic_energy
        instance.evolve_model(1.0 | nbody_system.time)
        energy_total_t1 = instance.potential_energy + instance.kinetic_energy
        
        self.assertAlmostRelativeEqual(energy_total_t0, energy_total_t1, 3)
        instance.stop()
        numpy.random.seed()
        
    
    def test17(self):
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
       
        instance = BHTree()
        instance.initialize_code()
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.collision_detection.enable()
        instance.evolve_model(0.1 | nbody_system.time)
        print instance.model_time
        self.assertTrue(instance.stopping_conditions.collision_detection.is_set())
        print instance.stopping_conditions.collision_detection.particles(0).key
        print instance.stopping_conditions.collision_detection.particles(1).key
        self.assertEquals(len(instance.stopping_conditions.collision_detection.particles(0)), 2 )
        p0 =  instance.stopping_conditions.collision_detection.particles(0)[0]
        p1 =  instance.stopping_conditions.collision_detection.particles(1)[0]
        self.assertNotEquals(p0, p1)
        print p0.x, p1.x
        self.assertTrue(p1.x - p0.x < 1.5| nbody_system.length)
       

        
    








