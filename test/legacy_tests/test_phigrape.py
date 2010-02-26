import os
import sys

from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.legacy.support import channel

from legacy_support import TestWithMPI
import path_to_test_results

import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestMPIInterface(TestWithMPI):
    
    def test0(self):
        instance = PhiGRAPEInterface()
        self.assertTrue("Refs not included yet" in instance.all_literature_references_string())

    def test1(self):
        instance = PhiGRAPEInterface()
        instance.setup_module()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_module()
        del instance
        
    def test2(self):
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.DDD
        instance = PhiGRAPEInterface()
        #channel.MessageChannel.DEBUGGER = None
        for x in [0.101, 4.0]:
            error = instance.set_eps2(x)
            self.assertEquals(error, 0)            
            value, error = instance.get_eps2()
            self.assertEquals(error, 0)
            self.assertEquals(x, value)
        del instance
        
    
    def test3(self):
        instance = PhiGRAPEInterface()
        instance.setup_module()
        
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
        instance.cleanup_module()
        del instance
    
    def test5(self):
        instance = PhiGRAPEInterface()
        instance.setup_module()
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
        instance.cleanup_module()
        
    def xtest6(self):
        instance = mpi_interface.PhiGRAPE()
        instance.setup_module()
        n = 4000
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
        for i in range(n-1):
            instance.add_particle(ids[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i]
                , values[i])
                
        retrieved_state = instance.get_state(1)
        self.assertEquals(1.0,  retrieved_state['mass'])
        instance.cleanup_module()

    def test7(self):
        instance = PhiGRAPEInterface()
        instance.setup_module()
        instance.set_eps2(0.0**2)
        instance.set_eta(0.01,0.02)

        instance.new_particle( 
            [1.0,1.0,1.0],
            [0.0,0.0,0.0],
            [1.0,0.0,-1.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0] )
        instance.initialize_particles(0.0) 
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']

        self.assertEqual( Ek, 0.5)
        self.assertEqual( Ep, -2.5)        
        instance.delete_particle(2)
        instance.reinitialize_particles() 
        n=instance.get_number_of_particles()['number_of_particles']
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']

        self.assertEqual( n, 2)
        self.assertEqual( Ek, 0.)
        self.assertEqual( Ep, -0.5)        

        instance.cleanup_module()
        del instance

    def test8(self):
        instance = PhiGRAPEInterface()
        instance.setup_module()
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
        instance.initialize_particles(0.0) 
        #HAS NO RESULT...
        result = instance.evolve(3.14159)  
        
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
        del instance
        

class TestSunAndEarthSystem(TestWithMPI):
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

        instance = PhiGRAPE(convert_nbody)
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.set_eta(0.01,0.02)
        instance.setup_module()
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]

        instance.setup_particles(stars)
        instance.initialize_particles(0.0)

        instance.evolve_model(365 | units.day)

        instance.update_particles(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.update_particles(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.update_particles(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        #self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_module()
        
        del instance
        

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = PhiGRAPE(convert_nbody)
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.set_eta(0.01,0.02)
        instance.setup_module()
        instance.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.setup_particles(stars)
        instance.initialize_particles(0.0)

        for x in range(1,365,1):
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
            output_file = os.path.join(test_results_path, "phiGRAPE-earth-sun2.svg")
            figure.savefig(output_file)
        
        instance.cleanup_module()
        del instance
        
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = PhiGRAPE(convert_nbody)
        #channel.MessageChannel.DEBUGGER = None

        instance.setup_module()
        
        particles = core.Particles(2)
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


        
