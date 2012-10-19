from amuse.units import units
from amuse.datamodel import Particles
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase

from amuse.ext.sticky_spheres import StickySpheres

class TestStickySpheres(TestCase):
    
    def test1(self):
        colliders = Particles(2)
        colliders.mass = [5, 2] | units.kg
        colliders.position = [[0.0, 0.0, 0.0], [0.7, 1.4, -0.35]] | units.m
        colliders.velocity = [[0.4, -0.6, 0.0], [0.0, 0.0, -3.0]] | units.m / units.s
        self.assertAlmostEqual(colliders.center_of_mass_velocity().length(), 1.0 | units.m / units.s)
        
        merged = StickySpheres().handle_collision(colliders[0], colliders[1])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(merged.mass, 7 | units.kg)
        self.assertAlmostEqual(merged.position, [0.2, 0.4, -0.1] | units.m)
        self.assertAlmostEqual(merged.velocity, ([2.0, -3.0, -6.0] | units.m / units.s) / 7.0)
        self.assertAlmostEqual(merged.velocity.length(), 1.0 | units.m / units.s)
        copy = colliders.copy()
        copy.move_to_center()
        self.assertAlmostEqual(colliders.kinetic_energy(), merged.as_set().kinetic_energy() + copy.kinetic_energy())
    
    def test2(self):
        colliders = Particles(2)
        colliders.mass = [5, 5] | units.kg
        colliders.position = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]] | units.m
        colliders.velocity = [[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]] | units.m / units.s
        
        merged = StickySpheres(mass_loss=0.2).handle_collision(colliders[0], colliders[1])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(merged.mass, 8 | units.kg)
        self.assertAlmostEqual(merged.position, [0.5, 0.5, 0.5] | units.m)
        self.assertAlmostEqual(merged.velocity, [1.0, 1.0, 1.0] | units.m / units.s)
        copy = colliders.copy()
        copy.move_to_center()
        self.assertAlmostEqual(colliders.kinetic_energy(), merged.as_set().kinetic_energy() / 0.8 + copy.kinetic_energy())
    
    def test3(self):
        self.assertRaises(AmuseException, StickySpheres, mass_loss=-0.1, expected_message=
            "Mass-loss fraction must be in the range [0, 1)")
        self.assertRaises(AmuseException, StickySpheres, mass_loss=1.0, expected_message=
            "Mass-loss fraction must be in the range [0, 1)")
        
        particles = Particles(6)
        particles.mass = range(1, 7) | units.kg
        particles.position = [[i, 1.0, 2.0] for i in range(1, 7)] | units.m
        particles.velocity = [[1.0, 0.0, 1.0], [0.0, -1.0, -1.0]] | units.m / units.s
        
        for fraction in [0.01, 0.1, 0.5]:
            sticky_spheres = StickySpheres(mass_loss=fraction)
            for i in range(0, 6, 2):
                colliders = particles[i:i+2]
                merged = sticky_spheres.handle_collision(colliders[0], colliders[1])
                self.assertTrue(isinstance(merged, Particles))
                self.assertAlmostEqual(merged.mass, 
                    (2 * i +3.0) * (1 - fraction) | units.kg)
                self.assertAlmostEqual(merged.position, 
                    [((i+1)**2 + (i+2)**2)/(2*i+3.0), 1.0, 2.0] | units.m)
                self.assertAlmostEqual(merged.velocity, 
                    ([i+1, -(i+2), -1.0] | units.m / units.s) / (2*i+3.0))
    

