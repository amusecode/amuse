import os
import sys
import time
from amuse.community.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE

from amuse.test.amusetest import TestWithMPI

import numpy
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.ic.plummer import new_plummer_sphere
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestMPIInterface(TestWithMPI):
    
    def test0(self):
        instance = PhiGRAPEInterface()
        self.assertTrue("Harfst, S., Gualandris, A., Merritt, D., Spurzem, R., Portegies Zwart, S., & Berczik, P." 
            in instance.all_literature_references_string())
        instance.stop()

    def test1(self):
        instance = PhiGRAPEInterface()
        instance.initialize_code()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        #instance.cleanup_code()
        instance.stop()
        
    def test2(self):
        instance = PhiGRAPEInterface()
        for x in [0.101, 4.0]:
            error = instance.set_eps2(x)
            self.assertEquals(error, 0)            
            value, error = instance.get_eps2()
            self.assertEquals(error, 0)
            self.assertEquals(x, value)
        instance.stop()
    
    def test3(self):
        instance = PhiGRAPEInterface()
        instance.initialize_code()
        
        instance.new_particle([11.0,12.0,13.0,14.0]
            , [2.0,3.0,4.0,5.0]
            , [2.1,3.1,4.1,5.1]
            , [2.2,3.2,4.2,5.2]
            , [2.3,3.3,4.3,5.3]
            , [2.4,3.4,4.4,5.4]
            , [2.5,3.5,4.5,5.5]
            , [2.6,3.6,4.6,5.6])
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        retrieved_state = instance.get_state([2,3,4])
        self.assertEquals(12.0,  retrieved_state['mass'][0])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 4)
        #instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        instance = PhiGRAPEInterface()
        instance.initialize_code()
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
        retrieved_state = instance.get_state(3999)
        self.assertEquals(3999.0,  retrieved_state['mass'])
        #instance.cleanup_code()
        instance.stop()

    def test6(self):
        instance = PhiGRAPEInterface()
        instance.initialize_code()
        n = 4000
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
        for i in range(n-1):
            instance.new_particle(
                  values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i])
                
        retrieved_state = instance.get_state(1)
        self.assertEquals(1.0,  retrieved_state['mass'])
        #instance.cleanup_code()
        instance.stop()

    def test7(self):
        instance = PhiGRAPEInterface()#(debugger="xterm")
        instance.initialize_code()
        
        instance.set_eps2(0.0**2)
        instance.set_eta(0.01,0.02)
        instance.commit_parameters()

        instance.new_particle( 
            [1.0,1.0,1.0],
            [0.0,0.0,0.0],
            [1.0,0.0,-1.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0] )
        instance.recommit_particles()
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        self.assertEqual( Ek, 0.5)
        self.assertEqual( Ep, -2.5)    
        instance.delete_particle(2)
        instance.recommit_particles()
        n=instance.get_number_of_particles()['number_of_particles']
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
    
        self.assertEqual( n, 2)
        self.assertEqual( Ek, 0.)
        self.assertEqual( Ep, -0.5)    
    
        #instance.cleanup_code()
        instance.stop()

    def xtest8(self):
        instance = PhiGRAPEInterface()
        instance.initialize_code()
        instance.set_eps2(0.0**2)
        instance.set_eta(0.01,0.02)
        instance.new_particle( 
            [0.01,0.01],
            [0.1,0.1],
            [10.,-10.],
            [0.0,0.0],
            [0.0,0.0],
            [-5.0,5.0],
            [0.0,0.0],
            [0.0,0.0] )
        instance.commit_particles() 
        #HAS NO RESULT...
        result = instance.evolve_model(3.14159)  
        
        tnow=instance.get_time()['time']
        print "after evolve(pi), tnow = %f" %  (tnow)
        #self.assertEqual( id1, 1)
        """
        instance.evolve(instance.get_time(),1)
        id2=instance.find_colliding_secondary(id1)
        self.assertEqual( id2, 2)
        self.assertAlmostEqual( tnow, 2.,2)
        state1 = instance.get_state(id1)
        state2 = instance.get_state(id2)
        self.assertTrue( abs(state1['x'] - state2['x'])<0.2)
        """
        instance.stop()
        
    
    def test8(self):
        instance = PhiGRAPEInterface()
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
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 10.0]) / 2.0)
        
    def test9(self):
        instance = PhiGRAPEInterface()
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
        

class TestPhigrape(TestWithMPI):
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
        instance = PhiGRAPE(convert_nbody)#, debugger="xterm")
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.set_eta(0.01,0.02)
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
    
        instance.particles.add_particles(stars)
        instance.commit_particles()
    
        position_at_start = earth.position.value_in(units.AU)[0]
        instance.evolve_model(365 | units.day)
    
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        #self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_code()
        
        instance.stop()

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = PhiGRAPE(convert_nbody)
        
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.set_eta(0.01,0.02)
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
            output_file = os.path.join(test_results_path, "phiGRAPE-earth-sun2.svg")
            figure.savefig(output_file)
        
        instance.cleanup_code()
        instance.stop()
        
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = PhiGRAPE(convert_nbody)

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
        instance = PhiGRAPE()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.set_eta(0.01,0.02)

        
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
        
    def test6(self):
        instance = PhiGRAPE()
        instance.initialize_code()
        
        particles = datamodel.Particles(6)
        particles.mass = nbody_system.mass.new_quantity(range(1,7))
        particles.radius =   0.00001 | nbody_system.length
        particles.position = [[-1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0],[0.0,0.0,1.0]] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        instance.commit_particles()
        copyof = instance.particles.copy()
        
        self.assertEquals(2 | nbody_system.mass, copyof[1].mass)  
        
        instance.stop()
        
    def test7(self):
        instance = PhiGRAPE()
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.set_eta(0.01,0.02)

        
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        instance.commit_particles()
        
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
        
    def test8(self):
        particles = datamodel.Particles(6)
        particles.mass = nbody_system.mass.new_quantity(range(1,7))
        particles.radius =   0.00001 | nbody_system.length
        particles.position = [[-1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0],[0.0,0.0,1.0]] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        
        for current_mode in ['g6lib','gpu','grape','pg']:
            try:
                instance = PhiGRAPE(mode = current_mode)
            except:
                print "Running PhiGRAPE with mode=", current_mode, " was unsuccessful."
            else:
                print "Running PhiGRAPE with mode=", current_mode, "... ",
                instance.initialize_code()
                instance.particles.add_particles(particles)
                instance.commit_particles()
                instance.evolve_model(0.1 | nbody_system.time)
                instance.cleanup_code()
                instance.stop()
                print "ok"
                
    
    def test9(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = PhiGRAPE(convert_nbody)#, debugger="xterm")
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.initialize_gpu_once = 1
        instance.set_eta(0.01,0.02)
        instance.dt_dia = 5000
        
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
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        #self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_code()
        
        instance.stop()

    def test10(self):
        instance = PhiGRAPE()
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        
        stars = new_plummer_sphere(100)
        stars.radius = 0 | nbody_system.length
        
        instance.particles.add_particles(stars)
        channel = stars.new_channel_to(instance.particles)
        
        instance.evolve_model(0.001 | nbody_system.time)
    
        e0 = instance.kinetic_energy + instance.potential_energy
        
        stars.mass *= 0.9
        channel.copy()
        
        instance.synchronize_model()
        
        e1 = instance.kinetic_energy + instance.potential_energy
        
        delta_e = e1 - e0
        
        self.assertTrue(e1 != e0)
        
        instance.stop()

    def test11(self):
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
       
        instance = PhiGRAPE()
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

    def test12(self):
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
       
        instance = PhiGRAPE()
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

    def test13(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass
       
        instance = PhiGRAPE()
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 2
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 2 | units.none)
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody_system.time)

        instance.stop()

    def test14(self):
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

        instance = PhiGRAPE()
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
        instance = PhiGRAPE(number_of_workers = 2)#, redirection = "none")
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.set_eta(0.01,0.02)
        instance.dt_dia = 5000
        print 1
        particles = datamodel.Particles(2)
        particles.x = [-1.0,1.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass
        
        instance.particles.add_particles(particles)
        instance.commit_particles()
    
        instance.evolve_model(0.01 | nbody_system.time)
    
        instance.particles.copy_values_of_all_attributes_to(particles)
        
        print instance.particles.position.x
        self.assertEquals(instance.particles[0].position.x, -instance.particles[1].position.x)
        self.assertAlmostRelativeEquals(instance.particles[1].position.x, 0.999969482111 | nbody_system.length, 8)
        
        instance.evolve_model(0.10 | nbody_system.time)
    
        instance.particles.copy_values_of_all_attributes_to(particles)
        
        print instance.particles.position.x
        self.assertEquals(instance.particles[0].position.x, -instance.particles[1].position.x)
        self.assertAlmostRelativeEquals(instance.particles[1].position.x, 0.99804560161 | nbody_system.length, 8)
        
        instance.evolve_model(0.50 | nbody_system.time)
    
        instance.particles.copy_values_of_all_attributes_to(particles)
        
        print instance.particles.position.x
        self.assertEquals(instance.particles[0].position.x, -instance.particles[1].position.x)
        self.assertAlmostRelativeEquals(instance.particles[1].position.x, 0.968416814302 | nbody_system.length, 8)
        
        instance.cleanup_code()
        
        instance.stop()
    
    
    def test16(self):
        instance = PhiGRAPE(number_of_workers = 2)
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.set_eta(0.01,0.02)
        instance.dt_dia = 5000
        print 1
        particles = datamodel.Particles(2)
        particles.x = [-1.0,1.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = [1.0, 0.5] | nbody_system.mass
        
        instance.particles.add_particles(particles)
        instance.commit_particles()
    
        instance.evolve_model(0.01 | nbody_system.time)
    
        instance.particles.copy_values_of_all_attributes_to(particles)
        
        print instance.particles.position.x
        self.assertAlmostRelativeEquals(instance.particles[0].position.x, -0.999984741095 | nbody_system.length, 8)
        self.assertAlmostRelativeEquals(instance.particles[1].position.x,  0.999969482189 | nbody_system.length, 8)
        
        instance.evolve_model(0.10 | nbody_system.time)
    
        instance.particles.copy_values_of_all_attributes_to(particles)
        
        print instance.particles.position.x
        self.assertAlmostRelativeEquals(instance.particles[0].position.x, -0.999022960148 | nbody_system.length, 8)
        self.assertAlmostRelativeEquals(instance.particles[1].position.x, 0.998045920305 | nbody_system.length, 8)
        
        instance.evolve_model(0.50 | nbody_system.time)
    
        instance.particles.copy_values_of_all_attributes_to(particles)
        
        print instance.particles.position.x
        self.assertAlmostRelativeEquals(instance.particles[0].position.x, -0.984250788742 | nbody_system.length, 8)
        self.assertAlmostRelativeEquals(instance.particles[1].position.x, 0.968501583383 | nbody_system.length, 8)
        
        instance.cleanup_code()
        
        instance.stop()
        
    def test17(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = PhiGRAPE(convert_nbody, number_of_workers=2)
        instance.initialize_code()
    
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.set_eta(0.01,0.02)
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
    
        instance.particles.add_particles(stars)
        instance.commit_particles()
    
        position_at_start = earth.position.value_in(units.AU)[0]
        print instance.particles[0].position.as_quantity_in(units.AU)
        print instance.particles[1].position.as_quantity_in(units.AU)
        instance.evolve_model(365 | units.day)
        print instance.particles[0].position.as_quantity_in(units.AU)
        print instance.particles[1].position.as_quantity_in(units.AU)
    
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 3)
        
        print instance.particles[0].position.as_quantity_in(units.AU)
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        print instance.particles[0].position.as_quantity_in(units.AU)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        #self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_code()
        
        instance.stop()

    
