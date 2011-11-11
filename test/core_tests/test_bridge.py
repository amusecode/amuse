from amuse.test import amusetest

from amuse.support.exceptions import AmuseException

import numpy
import time
import sys
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel
from amuse.couple import bridge

class TestCalculateFieldForParticles(amusetest.TestCase):
    
    def test1(self):
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        
        instance = bridge.CalculateFieldForParticles(particles = particles, gravity_constant = nbody_system.G)
        
        zero = 0.0 | nbody_system.length
        print instance.get_gravity_at_point([zero], [1.0] | nbody_system.length, [zero], [zero])
        fx, fy, fz = instance.get_gravity_at_point([zero], [1.0] | nbody_system.length, [zero], [zero])
        self.assertAlmostEqual(fx, [0.0] | nbody_system.acceleration, 6)
        self.assertAlmostEqual(fy, [0.0] | nbody_system.acceleration, 6)
        self.assertAlmostEqual(fz, [0.0] | nbody_system.acceleration, 6)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point([zero], [x0], [zero], [zero])
            potential1 = instance.get_potential_at_point([zero], [x1], [zero], [zero])
            fx0, fy0, fz0 = instance.get_gravity_at_point([zero], [x0], [zero], [zero])
            fx1, fy1, fz1 = instance.get_gravity_at_point([zero], [x1], [zero], [zero])
            
            self.assertAlmostEqual(fy0[0], 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fz0[0], 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fy1[0], 0.0 | nbody_system.acceleration, 6)
            self.assertAlmostEqual(fz1[0], 0.0 | nbody_system.acceleration, 6)
            
            self.assertAlmostEqual(fx0, -1.0 * fx1, 5)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0[0], 5)
            self.assertAlmostEqual(potential0, potential1, 6)
