import os
import sys
import numpy

from amuse.legacy.octgrav.interface import OctgravInterface, Octgrav

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.legacy import channel

from amuse.ext.plummer import *

from amuse.test.amusetest import TestWithMPI

class TestMPIInterface(TestWithMPI):

    def test1(self):
    
        instance = self.new_instance_of_an_optional_code(OctgravInterface)
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_code()
        instance.stop()

    def test2(self):
    
        instance = self.new_instance_of_an_optional_code(OctgravInterface)
        instance.initialize_code()
        instance.new_particle(
            [1,10,100],
            [2,11,101],
            [3,12,102],
            [4,13,103],
            [5,14,104],
            [6,15,105],
            [7,16,106],
            [8,17,107])
    
        particle1_state = instance.get_state(1)
        self.assertEquals(1,   particle1_state['mass'])
    
        particle2_state = instance.get_state(2)
        self.assertEquals(10,  particle2_state['mass'])
    
        instance.delete_particle(1)
    
        size_result = instance.get_number_of_particles()
        self.assertEquals(2, size_result['number_of_particles'])
    
        new_particle1_state = instance.get_state(1)
        self.assertEquals(10,  new_particle1_state['mass'])
    
        new_particle_result = instance.new_particle(
                                  1000,
                                  1001,
                                  1002,
                                  1003,
                                  1004,
                                  1005,
                                  1006,
                                  1007)
        self.assertEquals(4, new_particle_result['index_of_the_particle'],4)
    
        new_particle4_state = instance.get_state(4)
        self.assertEquals(1000,  new_particle4_state['mass'])
    
        instance.cleanup_code()
        instance.stop()

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

    def test0(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, units.AU)
        instance = self.new_instance_of_an_optional_code(Octgrav, convert_nbody)
        instance.initialize_code()
        self.assertAlmostRelativeEqual(0.01, instance.parameters.epsilon_squared.value_in(units.AU**2), 2)#default
        instance.parameters.epsilon_squared = 0.05 | units.AU**2
        self.assertAlmostRelativeEqual(0.05, instance.parameters.epsilon_squared.value_in(units.AU**2), 6)
       
        self.assertAlmostEqual(0.8|units.none, instance.parameters.opening_angle, 6)#default
        instance.parameters.opening_angle = 0.5
        self.assertAlmostEqual(0.5|units.none, instance.parameters.opening_angle, 6)
        instance.parameters.timestep = 1.0 |units.s
        self.assertEqual(1.0|units.s, instance.parameters.timestep)
        instance.stop()

    def test1(self):
        plummer_size = 500
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        plummer =  MakePlummerModel(plummer_size, convert_nbody)
        stars = plummer.result
        stars.radius = range(1, plummer_size+1)|units.km

        instance = self.new_instance_of_an_optional_code(Octgrav, convert_nbody)
        instance.particles.add_particles(stars)

        instance.evolve_model(5 | units.day)
        energy_total_init = instance.potential_energy + instance.kinetic_energy
        instance.evolve_model(100 | units.day)
        energy_total_final = instance.potential_energy + instance.kinetic_energy

        self.assertAlmostRelativeEqual(energy_total_init, energy_total_final, 2)
        instance.stop()

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        instance = self.new_instance_of_an_optional_code(Octgrav, convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        instance.parameters.stopping_conditions_number_of_steps = 1
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps,1|units.none)
    
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        instance.particles.add_particles(stars)
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(365.0 | units.day)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        instance.cleanup_code()
        instance.stop()
        
    def test3(self):
        particles = core.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass
       
        instance = self.new_instance_of_an_optional_code(Octgrav)
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 20
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 20|units.none)
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody_system.time)
        instance.stop()

    def test4(self):
        plummer_size = 500
        plummer =  MakePlummerModel(plummer_size)
        stars = plummer.result
        stars.radius=0|nbody_system.length

        instance = self.new_instance_of_an_optional_code(Octgrav)
        instance.particles.add_particles(stars)

        instance.synchronize_model()

        ax,ay,az=instance.get_gravity_at_point(0. | nbody_system.length,
                0. | nbody_system.length, 
                100. | nbody_system.length,
                0. | nbody_system.length)

        self.assertAlmostEqual(ax.number,0., 3)
        self.assertAlmostRelativeEqual(ay.number,-1./100**2, 3)
        self.assertAlmostEqual(az.number,0., 3)


        pot=instance.get_potential_at_point([0.,0.]|nbody_system.length,
                [0.,100] | nbody_system.length, 
                [100.,0.] | nbody_system.length,
                [0.,0.] | nbody_system.length)

        self.assertAlmostRelativeEqual(pot.number,[-1/100.,-1/100.], 3)

        instance.stop()
