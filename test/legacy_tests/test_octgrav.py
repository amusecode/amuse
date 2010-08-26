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
        instance.cleanup_module()
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

        instance.cleanup_module()
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
       
        self.assertAlmostEqual(0.8|units.none, instance.parameters.openings_angle, 6)#default
        instance.parameters.openings_angle = 0.5
        self.assertAlmostEqual(0.5|units.none, instance.parameters.openings_angle, 6)
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

        instance.evolve_model(1 | units.day)
        energy_total_init = instance.potential_energy + instance.kinetic_energy
        instance.evolve_model(100 | units.day)
        energy_total_final = instance.potential_energy + instance.kinetic_energy

        self.assertAlmostRelativeEqual(energy_total_init, energy_total_final, 2)
        instance.stop()

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Octgrav(convert_nbody)
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
        instance.cleanup_module()
        instance.stop()
        

