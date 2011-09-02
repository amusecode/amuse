import numpy
from amuse.test.amusetest import TestCase


from amuse.ext.solarsystem import new_solar_system, new_solar_system_for_mercury
from amuse.units import units
from amuse.units import constants
class TestSolarSystem(TestCase):
    
    def test1(self):
        print "Test 1: testing new_solar_system_for_mercury"
        sun, orbiters = new_solar_system_for_mercury()
        
        expected_attributes = set(["name", "mass", "radius", "j2", "j4", "j6", "Lx", "Ly", "Lz"])
        self.assertEqual(set(sun.get_attribute_names_defined_in_store()), expected_attributes)
        
        expected_attributes = set(["name", "mass", "radius", "density", "x", "y", "z", "vx", "vy", "vz", "Lx", "Ly", "Lz", "celimit"])
        self.assertEqual(set(orbiters.get_attribute_names_defined_in_store()), expected_attributes)
    
    def test2(self):
        print "Test 2: testing new_solar_system"
        particles = new_solar_system()
        print particles
        
        expected_attributes = set(["name", "mass", "radius", "x", "y", "z", "vx", "vy", "vz"])
        self.assertEqual(set(particles.get_attribute_names_defined_in_store()), expected_attributes)
        
        self.assertAlmostEqual(particles.center_of_mass(), [0, 0, 0] | units.m, in_units = units.AU)
        self.assertAlmostEqual(particles.center_of_mass_velocity(), [0, 0, 0] | units.m/units.s, in_units = units.AUd)
        # Particles are in center-of-mass(-velocity) coordinates, move them to heliocentric coordinates:
        self.assertTrue(particles[0].name == "SUN" | units.string)
        particles.position -= particles[0].position
        particles.velocity -= particles[0].velocity
        
        # Data from Carroll & Ostlie, An introduction to modern astrophysics, 1996
        eccentricity = [0.2056, 0.0068, 0.0167, 0.0934, 0.0483, 0.0560, 0.0461, 0.0097, 0.2482] | units.none
        semimajor_axis = [0.3871, 0.7233, 1.0000, 1.5237, 5.2028, 9.5388, 19.1914, 30.0611, 39.5294] | units.AU
        
        self.assertAlmostRelativeEqual(particles[2:-1].position.lengths_squared(), semimajor_axis[1:-1]**2, 1)
        
        # Somewhat more complicated test for more eccentric orbiters Mercury and Pluto:
        expected = (constants.G * (particles[1:].mass + particles[0].mass) * semimajor_axis * (1 - eccentricity**2)).sqrt()
        self.assertAlmostRelativeEqual(particles[1:].position.cross(particles[1:].velocity).lengths(), expected, 2)
    
