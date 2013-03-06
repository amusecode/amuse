import os
import sys
import time
from amuse.community.ph4.interface import ph4Interface, ph4

from amuse.test.amusetest import TestWithMPI

import numpy
import math
from amuse.units import nbody_system
from amuse.units import units, constants
from amuse import datamodel
from amuse.ic.plummer import new_plummer_model
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestMPIInterface(TestWithMPI):
    
    def test0(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.cleanup_code()
        instance.stop()
        
    def test1(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.set_eta(0.01)
        
        index, error = instance.new_particle(11.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(index, 1)
        error = instance.commit_particles()
        self.assertEquals(error, 0)
        
        retrieved_state = instance.get_state(index)
        self.assertEquals(retrieved_state['__result'], 0)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_code()
        instance.stop()
        
    def test2(self):
        instance = ph4Interface()
        instance.initialize_code()
        for x in [0.101, 4.0]:
            error = instance.set_eps2(x)
            self.assertEquals(error, 0)            
            value, error = instance.get_eps2()
            self.assertEquals(error, 0)
            self.assertEquals(x, value)
        instance.cleanup_code()
        instance.stop()
        
    
    def test3(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.set_eta(0.01)
        
        instance.new_particle([11.0,12.0,13.0,14.0]
            , [2.1,3.1,4.1,5.1]
            , [2.2,3.2,4.2,5.2]
            , [2.3,3.3,4.3,5.3]
            , [2.4,3.4,4.4,5.4]
            , [2.5,3.5,4.5,5.5]
            , [2.6,3.6,4.6,5.6]
            , [2.0,3.0,4.0,5.0])
        error = instance.commit_particles()
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        retrieved_state = instance.get_state([2,3,4])
        self.assertEquals(12.0,  retrieved_state['mass'][0])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 4)
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.set_eta(0.01)
        n = 4000
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
        instance.new_particle(
              values
            , values
            , values
            , values
            , values
            , values
            , values
            , values)
        error = instance.commit_particles()
        retrieved_state = instance.get_state(1)
        self.assertEquals(1.0,  retrieved_state['mass'])
        retrieved_state = instance.get_state(3999)
        instance.cleanup_code()
        instance.stop()
        self.assertEquals(3999.0,  retrieved_state['mass'])
        
    def test6(self):
        instance = ph4Interface()#(debugger="xterm")
        instance.initialize_code()
        instance.set_eps2(0.0**2)
        instance.set_eta(0.01)
        instance.commit_parameters()
        
        instance.new_particle( 
            [1.0,1.0,1.0],
            [1.0,0.0,-1.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0] )
        instance.commit_particles()
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        self.assertEqual( Ek, 0.5)
        self.assertEqual( Ep, -2.5)    
        instance.delete_particle(2)
        instance.recommit_particles()
        n=instance.get_number_of_particles()['number_of_particles']
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
    
        instance.cleanup_code()
        instance.stop()
    
        self.assertEqual( n, 2)
        self.assertEqual( Ek, 0.)
        self.assertEqual( Ep, -0.5)    

    def test7(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.set_eps2(0.1**2)
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2 + 0.1**2), 8)
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
    
        instance.cleanup_code()
        instance.stop()
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 10.0]) / 2.0)
        
    
    def test8(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.set_eps2(0)
        instance.commit_parameters()
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 1.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -1.0 / numpy.sqrt(2.0**2), 8)
        
        potential, errorcode = instance.get_potential(id2)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2), 8)
        
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
    
        instance.cleanup_code()
        instance.stop()
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 1.0]) / 2.0)
        
        
    def test9(self):
        instance = ph4Interface()
        instance.initialize_code()
        instance.set_eta(0.01)
        index, error = instance.new_particle(11.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 10)
        self.assertEquals(error, 0)
        self.assertEquals(index, 10)
        #index, error = instance.new_particle(12.0, 3.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10)
        #print index, error 
        # self.assertEquals(error, -1)
        #self.assertEquals(index, 10)
        error = instance.commit_particles()
        self.assertEquals(error, 0)
        
        retrieved_state = instance.get_state(index)
        self.assertEquals(retrieved_state['__result'], 0)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_code()
        instance.stop()

class TestPH4(TestWithMPI):
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
        
    
    def test0(self):
        instance = ph4()
        instance.stop()
        
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = ph4(convert_nbody)#, redirection="none")#, debugger="xterm")
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.timestep_parameter = 0.01
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
    
        instance.particles.add_particles(stars)
        instance.commit_particles()
    
        instance.evolve_model(365 | units.day)
    
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_code()
        
        instance.stop()

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = ph4(convert_nbody)
        
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.timestep_parameter = 0.01
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.particles.add_particles(stars)
        instance.commit_particles()
    
        for x in range(1,365,30):
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
            output_file = os.path.join(test_results_path, "ph4-earth-sun2.svg")
            figure.savefig(output_file)
        
        instance.cleanup_code()
        instance.stop()
        
    
    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = ph4(convert_nbody)

        instance.initialize_code()
        
        particles = datamodel.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        instance.commit_particles()
        
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        self.assertEquals(instance.get_mass(1), 17.0| units.kg) 
        self.assertEquals(instance.get_mass(2), 33.0| units.kg)  
        
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        instance = ph4()
        instance.initialize_code()
        
        particles = datamodel.Particles(6)
        particles.mass = nbody_system.mass.new_quantity(range(1,7))
        particles.radius =   0.00001 | nbody_system.length
        particles.position = [[-1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0],[0.0,0.0,1.0]] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        instance.commit_particles()
        copyof = instance.particles.copy()
        
        instance.cleanup_code()
        instance.stop()
        
        self.assertEquals(2 | nbody_system.mass, copyof[1].mass)  
        
        
    def test5(self):
        instance = ph4()
        instance.initialize_code()
        instance.parameters.manage_encounters = 2
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length ** 2
        particles = datamodel.Particles(6)
        particles.mass =  [0.01, 0.1,  0.1, 0.1, 0.1, 0.1] | nbody_system.mass
        particles.radius =   0.1 | nbody_system.length
        particles.position = [
            [-1.0,0.0,0.0], #first two close together
            [-1.2,0.0,0.0],
            [0.0,4.0,0.0], #rest far away
            [0.0,5.0,0.0],
            [0.0,6.0,0.0],
            [0.0,7.0,0.0]
        ] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        
        instance.particles.add_particles(particles)
        
        instance.commit_particles()
        instance.evolve_model(0.1  | nbody_system.time)
        
        instance.update_particle_set()
        
        self.assertEquals(len(instance.particles), 5)
        self.assertEquals(instance.particles.index_in_code, [3,4,5,6,10])
        self.assertEquals(instance.particles.mass, [0.1,0.1,0.1,0.1,0.11] | nbody_system.mass)
        
        self.assertEquals(len(particles), 6)
        instance.particles.synchronize_to(particles)
        self.assertEquals(len(particles), 5)
        self.assertEquals(particles.mass, [0.1,0.1,0.1,0.1,0.11] | nbody_system.mass)
        
        binary_energy1, error = instance.legacy_interface.get_binary_energy()
        self.assertEquals(error, 0)
        self.assertTrue(binary_energy1 < 0)
        
        binary_energy2 = instance.get_binary_energy()
        
        instance.cleanup_code()
        instance.stop()
        
        self.assertEquals(binary_energy2.value_in(nbody_system.energy), binary_energy1)
        
    def test6(self):
        instance = ph4()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2

        
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        instance.commit_particles()
        
        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, 1.0| nbody_system.length, zero, zero)
        self.assertAlmostEqual(fx, 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration, 3)
        for x in (0.25, 0.5, 0.75):
            x0 = x| nbody_system.length
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
      
        instance.stop()
        
    def xtest7(self):
        instance = ph4()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.set_eta(0.01,0.02)

        
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        instance.initialize_particles(0.0 )
        
        zero = [0.0, 0.0, 0.0] | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, [0.5, 1.0, 1.5] | nbody_system.length, zero, zero)
        self.assertAlmostRelativeEqual(fx[0], -3.55555555556 | nbody_system.acceleration, 5)
        self.assertAlmostRelativeEqual(fy[0], 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostRelativeEqual(fz[0], 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostRelativeEqual(fx[1], 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostRelativeEqual(fy[1], 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostRelativeEqual(fz[1], 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostRelativeEqual(fx[2], 3.55555555556 | nbody_system.acceleration, 5)
        self.assertAlmostRelativeEqual(fy[2], 0.0 | nbody_system.acceleration, 3)
        self.assertAlmostRelativeEqual(fz[2], 0.0 | nbody_system.acceleration, 3)
        
        n = 512
        x = nbody_system.length.new_quantity(numpy.linspace(0.1, 1.9, n))
        zero = nbody_system.length.new_quantity(numpy.zeros(n))
        fx, fy, fz = instance.get_gravity_at_point(zero, x, zero, zero)
        for i in range(n/2):
            self.assertAlmostRelativeEqual(fx[i], - fx[n - 1 - i], 5)
        
        instance.stop()
        
    def xtest8(self):
        particles = datamodel.Particles(6)
        particles.mass = nbody_system.mass.new_quantity(range(1,7))
        particles.radius =   0.00001 | nbody_system.length
        particles.position = [[-1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0],[0.0,0.0,1.0]] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        
        for current_mode in ['gpu']:
            try:
                instance = ph4(mode = current_mode)
            except:
                print "Running PhiGRAPE with mode=", current_mode, " was unsuccessful."
            else:
                print "Running PhiGRAPE with mode=", current_mode, "... ",
                instance.initialize_code()
                instance.particles.add_particles(particles)
                instance.initialize_particles(0.0)
                instance.evolve_model(0.1 | nbody_system.time)
                instance.cleanup_code()
                instance.stop()
                print "ok"
                
    

    def test10(self):
        instance = ph4()
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.parameters.timestep_parameter = 0.01
        
        stars = new_plummer_model(100)
        
        instance.particles.add_particles(stars)
        channel = stars.new_channel_to(instance.particles)
        
        instance.evolve_model(0.001 | nbody_system.time)
    
        e0 = instance.kinetic_energy + instance.potential_energy
        
        stars.mass *= 0.9
        channel.copy()
        
        instance.synchronize_model()
        
        e1 = instance.kinetic_energy + instance.potential_energy
        
        instance.cleanup_code()
        instance.stop()
        
        delta_e = e1 - e0
        
        self.assertTrue(e1 != e0)
    
    def test11(self):
        print "Testing PH4 collision_detection"
        particles = datamodel.Particles(7)
        particles.mass = 0.000001 | nbody_system.mass
        particles.radius = 0.01 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed
        
        instance = ph4()
        instance.initialize_code()
        instance.parameters.set_defaults()
        instance.parameters.epsilon_squared = 0 | nbody_system.length**2
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        
        for i in range(3): # PH4 can handle only one collision (=closest) at a time
            instance.evolve_model(1.0 | nbody_system.time)
            
            self.assertTrue(collisions.is_set())
            self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
            self.assertEquals(len(collisions.particles(0)), 1)
            self.assertEquals(len(collisions.particles(0)), 1)
            self.assertEquals(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 5 - i)
            self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius), True)
        
            sticky_merged = datamodel.Particles(len(collisions.particles(0)))
            sticky_merged.mass = collisions.particles(0).mass + collisions.particles(1).mass
            sticky_merged.radius = collisions.particles(0).radius
            for p1, p2, merged in zip(collisions.particles(0), collisions.particles(1), sticky_merged):
                merged.position = (p1 + p2).center_of_mass()
                merged.velocity = (p1 + p2).center_of_mass_velocity()
        
            print instance.model_time
            print instance.particles
            instance.particles.remove_particles(collisions.particles(0) + collisions.particles(1))
            instance.particles.add_particles(sticky_merged)
        
        instance.evolve_model(1.0 | nbody_system.time)
        print
        print instance.model_time
        print instance.particles
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 1.0 | nbody_system.time)
        self.assertEquals(len(collisions.particles(0)), 1)
        self.assertEquals(len(collisions.particles(1)), 1)
        self.assertEquals(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 2)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()

    def xtest11(self):
        particles = datamodel.Particles(2)
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
       
        instance = ph4()
        instance.initialize_code()
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.collision_detection.enable()
        instance.evolve_model(0.5 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.collision_detection.is_set())
        self.assertEquals(len(instance.stopping_conditions.collision_detection.particles(0)), 2 )
        p0 =  instance.stopping_conditions.collision_detection.particles(0)[0]
        p1 =  instance.stopping_conditions.collision_detection.particles(1)[0]
        self.assertNotEquals(p0, p1)
        self.assertTrue(p1.x - p0.x < 1.5| nbody_system.length)
        instance.stop()

    def xtest12(self):
        particles = datamodel.Particles(2)
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
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass
       
        instance = ph4()
        instance.initialize_code()
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.pair_detection.enable()
        instance.evolve_model(1.5 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.pair_detection.is_set())
        self.assertEquals(len(instance.stopping_conditions.pair_detection.particles(0)), 2 )
        p0 =  instance.stopping_conditions.pair_detection.particles(0)[0]
        p1 =  instance.stopping_conditions.pair_detection.particles(1)[0]
        self.assertNotEquals(p0, p1)
        self.assertTrue(p1.x - p0.x < 1.5| nbody_system.length)
        instance.stop()

    def xtest13(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass
       
        instance = ph4()
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 2
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 2)
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody_system.time)

        instance.stop()

    def xtest14(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        very_short_time_to_evolve = 1 | units.s
        very_long_time_to_evolve = 1e9 | nbody_system.time

        instance = ph4()
        instance.initialize_code()
        instance.parameters.stopping_conditions_timeout = very_short_time_to_evolve
        self.assertEquals(instance.parameters.stopping_conditions_timeout, very_short_time_to_evolve)
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.timeout_detection.enable()
        start = time.time()
        instance.evolve_model(very_long_time_to_evolve)
        end = time.time()
        self.assertTrue(instance.stopping_conditions.timeout_detection.is_set())
        self.assertTrue((end-start)<very_short_time_to_evolve.value_in(units.s) + 2)#2 = some overhead compensation

        instance.stop()
        
    
    def test15(self):
        instance = ph4()
        instance.initialize_code()
        
        
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        mass = instance.particles[0].mass
        
        instance.cleanup_code()
        instance.stop()
        
        self.assertEquals(mass, 1.0 | nbody_system.mass)
    
    def test16(self):
        print "Testing ph4 states"
        stars = new_plummer_model(100)
        black_hole = datamodel.Particle()
        black_hole.mass = 1.0 | nbody_system.mass
        black_hole.radius =  0.0 | nbody_system.length
        black_hole.position = [0.0, 0.0, 0.0] | nbody_system.length
        black_hole.velocity = [0.0, 0.0, 0.0] | nbody_system.speed
        
        print "First do everything manually:"
        instance = ph4()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.parameters.timestep_parameter = 0.01
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(stars)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.particles.remove_particle(stars[0])
        instance.particles.add_particle(black_hole)
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        instance.recommit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(0.001 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.synchronize_model()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print "initialize_code(), commit_parameters(), (re)commit_particles(), " \
            "synchronize_model(), and cleanup_code() should be called " \
            "automatically before editing parameters, new_particle(), get_xx(), and stop():"
        instance = ph4()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.parameters.timestep_parameter = 0.01
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(stars)
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.particles.remove_particle(stars[0])
        instance.particles.add_particle(black_hole)
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        mass = instance.particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(0.001 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        mass = instance.particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'STOPPED')
    


    def test17(self):
        print "Testing parameter defaults"
       
        instance = ph4()
        instance.parameters.epsilon_squared = 0.5  | nbody_system.length * nbody_system.length
        
        
        
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        self.assertAlmostRelativeEquals( instance.parameters.epsilon_squared ,  0.5  | nbody_system.length * nbody_system.length)
        self.assertAlmostRelativeEquals( instance.parameters.timestep_parameter ,  0.14)
        self.assertAlmostRelativeEquals( instance.parameters.use_gpu , 0)
        self.assertAlmostRelativeEquals( instance.parameters.manage_encounters , 4)
    
    def test18(self):
        print "Testing effect of ph4 parameter epsilon_squared"
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 3.0037e-6] | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * particles.total_mass() / (1.0 | units.AU)).sqrt()
        particles.rotate(0.0, 0.0, math.pi/4)
        particles.move_to_center()
        
        tan_initial_direction = particles[1].vy/particles[1].vx
        self.assertAlmostEquals(tan_initial_direction, math.tan(-math.pi/4))
        tan_final_direction =  []
        for log_eps2 in range(-5,6,2):
            instance = ph4(converter)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
            instance.parameters.timestep_parameter = 0.001
            instance.commit_parameters()
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            tan_final_direction.append(instance.particles[1].velocity[1]/
                instance.particles[1].velocity[0])
            instance.cleanup_code()
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEquals(tan_final_direction[0], math.tan(math.pi / 4.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEquals(tan_final_direction[-1], tan_initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(tan_final_direction[i+1]-tan_final_direction[i]) for i in range(len(tan_final_direction)-1)]
        self.assertEquals(delta[len(tan_final_direction)/2 -1], max(delta))
    
    def test19(self):
        print "Testing ph4 properties"
        particles = new_plummer_model(1000, do_scale=True)
        particles.position += [1, 2, 3] | nbody_system.length
        cluster_velocity = [4, 5, 6] | nbody_system.speed
        particles.velocity += cluster_velocity
        external_kinetic_energy = (0.5 | nbody_system.mass) * cluster_velocity.length_squared()
        
        instance = ph4()
        instance.particles.add_particles(particles)
        
        kinetic_energy = instance.kinetic_energy - external_kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 0.25 | nbody_system.energy, 10)
        self.assertAlmostRelativeEqual(potential_energy, -0.5 | nbody_system.energy, 10)
        self.assertAlmostRelativeEqual(instance.total_mass, 1.0 | nbody_system.mass, 10)
        self.assertAlmostRelativeEqual(instance.center_of_mass_position, 
            [1, 2, 3] | nbody_system.length, 10)
        self.assertAlmostRelativeEqual(instance.center_of_mass_velocity, 
            [4, 5, 6] | nbody_system.speed, 10)
        initial_total_energy = kinetic_energy + potential_energy
        
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostRelativeEqual(instance.model_time, 0.1 | nbody_system.time, 3)
        kinetic_energy = instance.kinetic_energy - external_kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy+potential_energy, -0.25 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(instance.total_mass, 1.0 | nbody_system.mass, 3)
        self.assertAlmostRelativeEqual(instance.center_of_mass_position, 
            [1.4, 2.5, 3.6] | nbody_system.length, 3)
        self.assertAlmostRelativeEqual(instance.center_of_mass_velocity, 
            [4, 5, 6] | nbody_system.speed, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test20(self):
        instance = ph4()
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
        potentials = [-4.57110767688 , -2.66662518639,-2.13331892165] | nbody_system.length**2 / (nbody_system.time**2) 
        for x, pot in zip((0.25, 0.5, 0.75), potentials):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            potential1 = instance.get_potential_at_point(zero, x1, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            fx1, fy1, fz1 = instance.get_gravity_at_point(zero, x1, zero, zero)
            
            self.assertAlmostRelativeEquals(fy0, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostRelativeEquals(fz0, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostRelativeEquals(fy1, 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostRelativeEquals(fz1, 0.0 | nbody_system.acceleration, 6)
            
            self.assertAlmostRelativeEquals(fx0, -1.0 * fx1, 5)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            print fx, fx0
            self.assertAlmostRelativeEquals(fx, fx0, 2)
            self.assertAlmostRelativeEquals(potential0, potential1, 5)
            self.assertAlmostRelativeEquals(potential0, pot, 5)
        instance.stop()
        
    def test21(self):
        particles = new_plummer_model(200)
        particles.scale_to_standard()
        instance = ph4()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000 | nbody_system.length**2
        instance.particles.add_particles(particles)
        
        x = numpy.arange(-1,1,0.1) | nbody_system.length
        zero = numpy.zeros(len(x)) | nbody_system.length
        potential0 =  instance.get_potential_at_point(zero, x , zero, zero)
        instance.stop()
        for n in (2, 3, 4):
            
            instance = ph4(number_of_workers = n)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 0.00000 | nbody_system.length**2
            instance.particles.add_particles(particles)
            potential =  instance.get_potential_at_point(zero, x , zero, zero)
        
            self.assertAlmostRelativeEquals(potential0, potential, 8)
            instance.stop()
    
    def test22(self):
        particles = new_plummer_model(200)
        particles.scale_to_standard()
        try:
            instance = ph4(mode="gpu")
        except:
            self.skip("gpu mode is not available")
            
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000 | nbody_system.length**2
        instance.particles.add_particles(particles)
        
        x = numpy.arange(-1,1,0.1) | nbody_system.length
        zero = numpy.zeros(len(x)) | nbody_system.length
        gpu_potential =  instance.get_potential_at_point(zero, x , zero, zero)
        instance.stop()
        
        instance = ph4()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000 | nbody_system.length**2
        instance.particles.add_particles(particles)
        nogpu_potential =  instance.get_potential_at_point(zero, x , zero, zero)
    
        self.assertAlmostRelativeEquals(nogpu_potential, gpu_potential, 5)
        instance.stop()
    
    def test23(self):
        
        particles = datamodel.Particles(
            mass=[1,2] | nbody_system.mass,
            x=[-1,1] | nbody_system.length,
            y=[-1,1] | nbody_system.length,
            z=[-1,1] | nbody_system.length,
            vx=[-1,1] | nbody_system.speed,
            vy=[-1,1] | nbody_system.speed,
            vz=[-1,1] | nbody_system.speed
        )
        
        
        instance=ph4()
        
        overlay = datamodel.ParticlesOverlay(instance.particles)
        
        overlay.add_particles(particles)
        all_attributes = overlay.get_values_in_store(overlay.get_all_indices_in_store(), ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz'])
        
        self.assertEquals(all_attributes[0], [1,2] | nbody_system.mass )
        self.assertEquals(instance.particles.mass, [1,2] | nbody_system.mass )
        self.assertEquals(overlay.mass, [1,2] | nbody_system.mass )
        self.assertEquals(overlay.position, [[-1., -1., -1.], [ 1.,  1.,  1.]]  | nbody_system.length )

   
    def test24(self):
        
        particles = datamodel.Particles(
            mass=[1,2] | nbody_system.mass,
            x=[-1,1] | nbody_system.length,
            y=[-1,1] | nbody_system.length,
            z=[-1,1] | nbody_system.length,
            vx=[-1,1] | nbody_system.speed,
            vy=[-1,1] | nbody_system.speed,
            vz=[-1,1] | nbody_system.speed
        )
        
        instance=ph4()
        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        self.assertAlmostRelativeEquals(instance.particles.mass, [1,2] | nbody_system.mass)
        
        instance.cleanup_code()
        instance.initialize_code()
        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
        self.assertAlmostRelativeEquals(instance.particles.mass, [1,2] | nbody_system.mass)


    def test25(self):
        
        particles = datamodel.Particles(
            mass=[1,2] | nbody_system.mass,
            x=[-1,1] | nbody_system.length,
            y=[-1,1] | nbody_system.length,
            z=[-1,1] | nbody_system.length,
            vx=[-1,1] | nbody_system.speed,
            vy=[-1,1] | nbody_system.speed,
            vz=[-1,1] | nbody_system.speed,
            id=[3,4]
        )
        
        instance=ph4()
        
        instance.particles.add_particles(particles)
        self.assertEquals(instance.particles.index_in_code, [3,4])
        instance.commit_particles()
        instance.particles.remove_particle(particles[1])
        instance.recommit_particles()
        self.assertEquals(instance.particles.index_in_code, [3])
        instance.particles.add_particle(particles[1])
        instance.recommit_particles()
        self.assertEquals(instance.particles.index_in_code, [3,4])
           

    def test26(self):
        
        particles = datamodel.Particles(
            mass=[1,2] | nbody_system.mass,
            x=[-1,1] | nbody_system.length,
            y=[ 0,0] | nbody_system.length,
            z=[ 0,0] | nbody_system.length,
            vx=[-1,1] | nbody_system.speed,
            vy=[-1,1] | nbody_system.speed,
            vz=[-1,1] | nbody_system.speed,
            id=[3,4]
        )
        
        instance=ph4()
        
        instance.particles.add_particles(particles)
        self.assertEquals(instance.particles.index_in_code, [3,4])
        instance.commit_particles()
        print instance.particles[0].potential_in_code
        self.assertEquals(instance.particles[0].potential(gravitationalConstant=nbody_system.G), -(nbody_system.G * (2 | nbody_system.mass)) / (2 | nbody_system.length))
        self.assertEquals(instance.particles[0].potential_in_code, -(nbody_system.G * (2 | nbody_system.mass)) / (2 | nbody_system.length))
