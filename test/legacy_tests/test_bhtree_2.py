# -*- coding: utf-8 -*-
from legacy_support import TestWithMPI
import os
import sys

from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from support import path_to_test_results

import numpy

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
        instance.setup_module()

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

        self.assertEquals(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        self.assertEquals(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle']) 
        
        instance.cleanup_module()
        del instance

    def test2(self):
        instance = BHTreeInterface()
        instance.setup_module()

        for i in [1, 2, 3]:
            temp_particle = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEquals(i, temp_particle['index_of_the_particle'])

        self.assertEquals(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        self.assertEquals(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle']) 
        self.assertEquals(3, instance.get_index_of_next_particle(2)['index_of_the_next_particle'])
            
        instance.delete_particle(1)
      
        self.assertEquals(2, instance.get_number_of_particles()['number_of_particles'])
        
        #the deletion does a swap, so 3 is copied to 1, (overwriting old 1 and treesize -> treesize-1
        self.assertEquals(3, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEquals(1, instance.get_index_of_next_particle(2)['__result'])

        instance.cleanup_module()
        del instance
        
    def test3(self):
        interface = BHTreeInterface()
        interface.eps2 = 0.101
        self.assertEquals(0.101, interface.eps2)
        interface.eps2 = 0.110
        self.assertEquals(0.110, interface.eps2)
        del interface

    def test4(self):
        interface = BHTreeInterface()
        interface.flag_collision = 1
        self.assertEquals(1, interface.flag_collision)
        interface.flag_collision = 0
        self.assertEquals(0, interface.flag_collision)
        del interface

    def test5(self):
        interface = BHTreeInterface()
        interface.setup_module()
        
        interface.new_particle([10,20],[1,1],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0])
        retrieved_state = interface.get_state(1)
        
        self.assertEquals(10.0,  retrieved_state['mass'])
        self.assertEquals(1, retrieved_state['radius'])

        retrieved_state = interface.get_state([1,2])
        self.assertEquals(20.0,  retrieved_state['mass'][1])
        self.assertEquals(interface.get_number_of_particles()['number_of_particles'], 2)
        interface.cleanup_module()
        del interface

    def test6(self):
        interface = BHTreeInterface()
        interface.setup_module()
        interface.cleanup_module()
        del interface

class TestAmuseInterface(TestWithMPI):
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = BHTree(convert_nbody)
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        instance.setup_module()
        
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

        instance.setup_particles(stars)

        instance.evolve_model(365.0 | units.day)
        instance.update_particles(stars)
        
        postion_at_start = earth.position.value_in(units.AU)[0]
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
        instance.cleanup_module()
        del instance
        
    def test2(self):
	#not completed 
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = BHTree(convert_nbody)
        #instance.dt_dia = 1
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        #instance.timestep = 0.0001
        #instance.use_self_gravity = 0
        instance.setup_module()
        
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

        instance.setup_particles(stars)
    
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
               
            
            test_results_path = path_to_test_results.get_path_to_test_results()
            output_file = os.path.join(test_results_path, "bhtree-earth-sun.svg")
            figure.savefig(output_file)    
        
        instance.cleanup_module()
        del instance

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = BHTree(convert_nbody)
        #instance.dt_dia = 1
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        #instance.timestep = 0.0001
        #instance.use_self_gravity = 0
        instance.setup_module()
        
        
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
        
        instance.setup_particles(stars)
    
        for x in range(1,200,1):
            instance.evolve_model(x | units.day)
            instance.update_particles(stars)
            print instance.get_indices_of_colliding_particles()
            #print stars[0].position-stars[1].position
            stars.savepoint()
        
         

