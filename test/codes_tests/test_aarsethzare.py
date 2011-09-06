from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.aarsethzare.interface import AarsethZareInterface, AarsethZare

from amuse.ext import plummer
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class TestAarsethZareInterface(TestWithMPI):

    def test0(self):
        instance = AarsethZareInterface()
        instance.stop()
        
    def test1(self):
        
        instance = AarsethZareInterface()
        time, x, y, z, vx, vy, vz, error = instance.evolve_triple(
            [0,0,0],
            [1,1,1],
            [-10, 0, 10],
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [1.0,1.0,1.0]
        )
        
        print time, x, y, z, vx, vy, vz, error
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(y, [0,0,0])
        self.assertAlmostRelativeEquals(z, [0,0,0])
        self.assertAlmostRelativeEquals(x[-1], -x[0], 10)
        self.assertAlmostRelativeEquals(x[1], 0, 10)
        self.assertAlmostRelativeEquals(vx[-1], -vx[0], 10)
        self.assertAlmostRelativeEquals(vx[1], 0, 10)
        self.assertAlmostRelativeEquals(vy, [0,0,0])
        self.assertAlmostRelativeEquals(vz, [0,0,0])
        
    
    def test2(self):
        
        instance = AarsethZareInterface()
        time, x, y, z, vx, vy, vz, error = instance.evolve_triple(
            [0,0],
            [1,1],
            [-10, 0],
            [0,0],
            [0,0],
            [0,0],
            [0,0],
            [0,0],
            [1.0]
        )
        
        self.assertEquals(error, -1)

class TestAarsethZare(TestWithMPI):
    def new_system_of_sun_and_earth_and_moon(self):
        stars = datamodel.Stars(3)
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
        
        moon = stars[2]
        moon.mass = units.kg(7.3477e22 )
        moon.radius = units.km(1737.10) 
        moon.position = units.km(numpy.array((149.5e6 + 384399.0 ,0.0,0.0)))
        moon.velocity = ([0.0,1.022,0] | units.km/units.s) + earth.velocity
        return stars
        
    def test0(self):
        instance = AarsethZare()
        instance.stop()
    
    def test1(self):
        instance = AarsethZare()
        time, x, y, z, vx, vy, vz = instance.evolve_triple(
            [0,0,0] | nbody_system.time,
            [1,1,1] | nbody_system.mass,
            [-10, 0, 10] | nbody_system.length,
            [0,0,0] | nbody_system.length,
            [0,0,0] | nbody_system.length,
            [0,0,0] | nbody_system.speed,
            [0,0,0] | nbody_system.speed,
            [0,0,0] | nbody_system.speed,
            [1.0,1.0,1.0]| nbody_system.time
        )
        
        self.assertAlmostRelativeEquals(y, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(z, [0,0,0] | nbody_system.length )
        self.assertAlmostRelativeEquals(x[-1], -x[0], 10 )
        self.assertAlmostRelativeEquals(x[1], 0 | nbody_system.length , 10)
        self.assertAlmostRelativeEquals(vx[-1], -vx[0], 10)
        self.assertAlmostRelativeEquals(vx[1], 0 | nbody_system.speed, 10)
        self.assertAlmostRelativeEquals(vy, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(vz, [0,0,0] | nbody_system.speed)
        
    def test2(self):
        instance = AarsethZare()
        
        particles = datamodel.Particles(3)
        particles.mass = 1.0 | nbody_system.mass
        particles.position = [
                [-10.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
            ] | nbody_system.length
        particles.velocity = [0.0, 0.0, 0.0] | nbody_system.speed
        
        instance.particles.add_particles(particles)
        instance.evolve_model(1.0 | nbody_system.time)
        
        
        self.assertAlmostRelativeEquals(instance.particles.y, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.particles.z, [0,0,0] | nbody_system.length )
        self.assertAlmostRelativeEquals(instance.particles.x[-1], -instance.particles.x[0], 10 )
        self.assertAlmostRelativeEquals(instance.particles.x[1], 0 | nbody_system.length , 10)
        self.assertAlmostRelativeEquals(instance.particles.vx[-1], -instance.particles.vx[0], 10)
        self.assertAlmostRelativeEquals(instance.particles.vx[1], 0 | nbody_system.speed, 10)
        self.assertAlmostRelativeEquals(instance.particles.vy, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(instance.particles.vz, [0,0,0] | nbody_system.speed)
            
    
    
    def test3(self):
        instance = AarsethZare(nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.MSun))
        
        particles = self.new_system_of_sun_and_earth_and_moon()
        instance.particles.add_particles(particles)
        
        earth = instance.particles[1]
        
        position_at_start = earth.position[0]
        
        instance.evolve_model(365.0 | units.day)
        
        position_after_full_rotation = earth.position[0]
        
        self.assertAlmostRelativeEquals(position_at_start, position_after_full_rotation, 4)

        instance.evolve_model(1.5 | units.yr)
        
        position_after_half_a_rotation = earth.position[0]
        
        self.assertAlmostRelativeEquals(-position_at_start, position_after_half_a_rotation, 2)
        
        instance.evolve_model(1.75 | units.yr)
        
        position_after_quarter_a_rotation = earth.position[1]
        
        self.assertAlmostRelativeEquals(position_at_start, position_after_quarter_a_rotation, 2)
    
        instance.stop()
    
