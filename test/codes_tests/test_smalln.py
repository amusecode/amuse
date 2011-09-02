import os
import sys
import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.community.smallN.muse_dynamics_mpi import SmallNInterface, SmallN



from amuse.support.data import core

from mpi4py import MPI
from amuse.units import nbody_system
from amuse.units import units
class TestSmallNInterface(TestWithMPI):
    
    def test0(self):
        instance=self.new_instance(SmallNInterface)
        instance.initialize_code()
        instance.stop()

    def test1(self):
        instance=self.new_instance(SmallNInterface)
        instance.initialize_code()

        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEquals(0, res1)
        self.assertEquals(1, res2)

        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)
        print retrieved_state1
        self.assertEquals(11.0,  retrieved_state1[0])
        self.assertEquals(21.0,  retrieved_state2[0])
        self.assertEquals(0.0,  retrieved_state1[1])
        self.assertEquals(10.0,  retrieved_state2[1])

        instance.cleanup_code()
        instance.stop()
        
    
class TestSmallN(TestWithMPI):
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
        instance=self.new_instance(SmallN, convert_nbody)
        instance.initialize_code()
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
             
        instance.particles.add_particles(stars)
        delta_position = instance.particles[1].position - instance.particles[0].position
        r0 = delta_position.length()
        instance.evolve_model()
        
        delta_position_2 = instance.particles[1].position - instance.particles[0].position
        r1 = delta_position_2.length()
        
        self.assertAlmostRelativeEquals(r0, r1, 8)
        
        instance.cleanup_code()
        instance.stop()
