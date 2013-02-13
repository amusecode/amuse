import numpy
from amuse.test.amusetest import TestCase
from amuse.units import units
from amuse.datamodel import Particles
from amuse.ext.particles_with_color import *


class TestParticlesWithColor(TestCase):
    
    def test1(self):
        print "Test new_particles_with_color"
        with_color = new_particles_with_color(
            Particles(mass=[0.1, 1, 2, 3]|units.MSun), 
            lambda mass: numpy.select([mass < 2.5|units.MSun], [1]),
            lambda mass: numpy.select([numpy.logical_and(0.5|units.MSun < mass, mass < 2.5|units.MSun)], [1]),
            lambda mass: numpy.select([mass > 1.5|units.MSun], [1])
        )
        self.assertEqual(with_color.red,   [1, 1, 1, 0])
        self.assertEqual(with_color.green, [0, 1, 1, 0])
        self.assertEqual(with_color.blue,  [0, 0, 1, 1])
        self.assertEqual(with_color.color, [[1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 0, 1]])
        with_color.mass = [3, 2, 1, 0.1] | units.MSun
        self.assertEqual(with_color.red,   [0, 1, 1, 1])
        self.assertEqual(with_color.green, [0, 1, 1, 0])
        self.assertEqual(with_color.blue,  [1, 1, 0, 0])
        self.assertEqual(with_color.color, [[0, 0, 1], [1, 1, 1], [1, 1, 0], [1, 0, 0]])
    
    def test2(self):
        print "Test new_particles_with_blackbody_color from temperature"
        original = Particles(temperature=[100, 1000, 2000, 5000, 10000, 40000]|units.K)
        with_color = new_particles_with_blackbody_color(original)
        self.assertAlmostEqual(with_color.red,   [1.0000, 1.0000, 1.0000, 1.0000, 0.6324, 0.3565], 3)
        self.assertAlmostEqual(with_color.green, [0.0337, 0.0337, 0.2434, 0.7872, 0.7081, 0.4747], 3)
        self.assertAlmostEqual(with_color.blue,  [0.0000, 0.0000, 0.0000, 0.5797, 1.0000, 1.0000], 3)
        self.assertEqual(with_color.color.shape, (6, 3))
    
    def test3(self):
        print "Test new_particles_with_blackbody_color from internal energy"
        original = Particles(u=u_from_T([100, 1000, 2000, 5000, 10000, 40000]|units.K))
        with_color = new_particles_with_blackbody_color(original)
        self.assertAlmostEqual(with_color.red,   [1.0000, 1.0000, 1.0000, 1.0000, 0.6324, 0.3565], 3)
        self.assertAlmostEqual(with_color.green, [0.0337, 0.0337, 0.2434, 0.7872, 0.7081, 0.4747], 3)
        self.assertAlmostEqual(with_color.blue,  [0.0000, 0.0000, 0.0000, 0.5797, 1.0000, 1.0000], 3)
        self.assertEqual(with_color.color.shape, (6, 3))
    

