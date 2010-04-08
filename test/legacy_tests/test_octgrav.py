import os
import sys

from amuse.legacy.octgrav.interface import OctgravInterface, Octgrav

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

    def test1(self):
        #channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = self.new_instance_of_an_optional_code(OctgravInterface)
        if instance is None:
            return

        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_module()
        del instance

    def test2(self):
        instance = self.new_instance_of_an_optional_code(OctgravInterface)
        if instance is None:
            return
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
        del instance

class TestAmuseInterface(TestWithMPI):

    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        channel.MessageChannel.DEBUGGER = channel.MessageChannel.XTERM
        instance = Octgrav(convert_nbody)
        instance.parameters.epsilon_squared = 0.001 | units.AU**2
        #instance.parameters.openings_angle = 0.1 |units.none
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
        
        instance.particles.add_particles(stars)
        
        postion_at_start = earth.position.value_in(units.AU)[0]
        print postion_at_start
        instance.evolve_model(365.0 | units.day)
        instance.update_particles(stars)
        
        postion_after_full_rotation = earth.position.value_in(units.AU)[0]

        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 3)
        print postion_after_full_rotation
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

