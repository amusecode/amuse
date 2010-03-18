from legacy_support import TestWithMPI
import os
import sys

from amuse.legacy.hermite0.interface import HermiteInterface, Hermite

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
import path_to_test_results

import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestMPIInterface(TestWithMPI):

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

        instance.cleanup_module()
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
        del hermite

    def test4(self):
        hermite = HermiteInterface()
        hermite.flag_collision = 1
        self.assertEquals(1, hermite.flag_collision)
        hermite.flag_collision = 0
        self.assertEquals(0, hermite.flag_collision)
        del hermite

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
        hermite.cleanup_module() 

    def test6(self):
        hermite = HermiteInterface()
        hermite.initialize_code()
        
        hermite.new_particle([10,10],[1,1],[-1,1],[0,0], [0,0], [0,0], [0,0], [0,0])
        retrieved_state = hermite.get_state(0)
        
        retr = hermite.get_potential_at_point(0.01, 0,0,0)
        self.assertEqual(retr['phi'], -20.0)
        hermite.cleanup_module()
       

class TestAmuseInterface(TestWithMPI):
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
        hermite.parameters.epsilon_squared = 0.0 | units.AU**2
        hermite.initialize_code()
        hermite.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        hermite.setup_particles(stars)
        
        hermite.evolve_model(365.0 | units.day)
        hermite.update_particles(stars)
        print stars[0].x
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        hermite.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        hermite.update_particles(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        hermite.update_particles(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        hermite.cleanup_module()
        
        del hermite
        

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Hermite(convert_nbody)
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.initialize_code()
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
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
            output_file = os.path.join(test_results_path, "hermite-earth-sun2.svg")
            figure.savefig(output_file)
        
        
        
        instance.cleanup_module()
        del instance

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Hermite(convert_nbody)
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.initialize_code()
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
        
        instance.setup_particles(stars)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.update_particles(stars)
            stars.savepoint()
        
    
    def test4(self):
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


         
