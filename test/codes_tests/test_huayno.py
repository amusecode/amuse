from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.huayno.interface import HuaynoInterface, Huayno

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.ic import plummer
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestHuaynoInterface(TestWithMPI):
    
    def test1(self):
        instance = HuaynoInterface()
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
        instance = HuaynoInterface()
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
        huayno = HuaynoInterface()
        huayno.eps2 = 0.101
        self.assertEquals(0.101, huayno.eps2)
        huayno.eps2 = 0.110
        self.assertEquals(0.110, huayno.eps2)
        huayno.stop()

    def test5(self):
        huayno = HuaynoInterface()
        huayno.initialize_code()
        
        huayno.new_particle([10,20],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0],[1,1])
        retrieved_state = huayno.get_state(0)
        
        self.assertEquals(10.0,  retrieved_state['mass'])
        self.assertEquals(1, retrieved_state['radius'])
    
        retrieved_state = huayno.get_state([0,1])
        self.assertEquals(20.0,  retrieved_state['mass'][1])
        self.assertEquals(huayno.get_number_of_particles()['number_of_particles'], 2)
        huayno.cleanup_code() 
       

class TestHuayno(TestWithMPI):
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
    
        huayno = Huayno(convert_nbody)
        huayno.initialize_code()
        huayno.parameters.epsilon_squared = 0.0 | units.AU**2
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        
        huayno.particles.add_particles(stars)
        
        huayno.evolve_model(365.0 | units.day)
        huayno.particles.copy_values_of_all_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        huayno.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        huayno.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        huayno.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        huayno.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        huayno.cleanup_code()
        
        huayno.stop()
        

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        instance = Huayno(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        
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
            output_file = os.path.join(test_results_path, "huayno-earth-sun2.svg")
            figure.savefig(output_file)
              
        instance.cleanup_code()
        del instance

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Huayno(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | units.AU**2
        
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
        
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Huayno(convert_nbody)
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
        
        
        self.assertEquals(instance.get_mass(0), 17.0| units.kg) 
        self.assertEquals(instance.get_mass(1), 33.0| units.kg)  

    def test5(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Huayno(convert_nbody)
        instance.initialize_code()
        
        particles = datamodel.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        
        instance.set_state(1, 16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms)
        
        curr_state =  instance.get_state(1)
        for expected, actural in zip((16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms, 0 | units.m), curr_state):
            self.assertAlmostRelativeEquals(actural,expected)
        
        instance.set_state(1, 16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms , 20.0|units.m)
        
        curr_state =  instance.get_state(1)
        for expected, actural in zip((16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms, 20 | units.m), curr_state):
            self.assertAlmostRelativeEquals(actural,expected)
        
        
    
    def test6(self):
        print "Test6: Testing Huayno parameters"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = Huayno(convert_nbody)
        
        (value, error) = instance.legacy_interface.get_eps2()
        self.assertEquals(0, error)
        self.assertEquals(0.0, value)
        self.assertAlmostEquals(0.0 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEquals(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
                
        (value, error) = instance.legacy_interface.get_time()
        self.assertEquals(0, error)
        self.assertEquals(0.0, value)
        self.assertAlmostEquals(0.0 | units.yr, instance.parameters.begin_time, in_units=units.yr)
        for x in [1.0, 10.0, 100.0]:
            instance.parameters.begin_time = x | units.yr
            self.assertAlmostEquals(x | units.yr, instance.parameters.begin_time, in_units=units.yr)
        instance.stop()
    
    def test7(self):
        print "Test7: Testing effect of Huayno parameter epsilon_squared"
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
            instance = Huayno(convert_nbody)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            final_direction.append(math.atan((instance.particles[1].velocity[0]/
                instance.particles[1].velocity[1])))
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEquals(abs(final_direction[0]), abs(initial_direction+math.pi/2.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEquals(final_direction[-1], initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(final_direction[i+1]-final_direction[i]) for i in range(len(final_direction)-1)]
        self.assertEquals(delta[len(final_direction)/2 -1], max(delta))
        
        
    def test13(self):
        particles = plummer.new_plummer_model(31)
       
        instance = Huayno(number_of_workers=1)#, debugger="xterm")
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.01 | nbody_system.length ** 2
        instance.particles.add_particles(particles)
        
        instance.evolve_model(0.1 | nbody_system.time)
        instance.synchronize_model()
        expected_positions = instance.particles.position
        instance.stop()
    
    def test14(self):
        import hashlib

        numpy.random.seed(123456)
        particles = plummer.new_plummer_model(64)
        sha=hashlib.sha1()

        for itype in sorted(Huayno.inttypes._list()):
            if itype in ("KEPLER"): continue
            instance = Huayno()
            instance.parameters.inttype_parameter=getattr(Huayno.inttypes,itype)
            instance.particles.add_particles(particles)
            instance.evolve_model(0.125 | nbody_system.time)
            part_out= instance.particles.copy()
            sha.update(part_out.position.number.data)
            instance.stop()
        
        # this result is probably dependent on system architecture hence no good for assert
        print sha.hexdigest()
        print "9714521156a4d4befc1e414f75aa650e1f8836ae"
    
