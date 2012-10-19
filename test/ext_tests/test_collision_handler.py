from amuse.units import units
from amuse.datamodel import Particles, Particle
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase

from amuse.ext.sticky_spheres import StickySpheres
from amuse.couple.collision_handler import CollisionHandler


class CollisionCodeForTesting(object):
    
    class ParametersForTesting(object):
        def __init__(self):
            self.mass_unit = units.kg
    
    stellar_evolution_code_required = False
    gravity_code_required = False
    
    def __init__(self, next_mass = 1 | units.kg):
        self.next_mass = next_mass
        self.parameters = self.ParametersForTesting()
    
    def handle_collision(self, primary, secondary, stellar_evolution_code=None):
        result = Particles(1)
        result.mass = self.next_mass.as_quantity_in(self.parameters.mass_unit)
        self.next_mass += 1 | units.kg
        
        if not stellar_evolution_code is None:
            se_colliders = (primary + secondary).get_intersecting_subset_in(stellar_evolution_code.particles)
            result.radius = se_colliders.radius.sum()
            
            def internal_structure(set, particle=None):
                return dict(mass=result.mass, radius=result.radius)
            
            result.add_function_attribute("internal_structure", None, internal_structure)
        
        return result
    

class GravityCodeForTesting(object):
    
    def __init__(self):
        self.particles = Particles(6)
        self.particles.mass = 1 | units.MSun
        self.particles.radius = 3 | units.RSun
        self.particles.position = [[i,2*i,3*i] for i in range(6)] | units.AU
        self.particles.velocity = [[i,i**2,i**3] for i in range(6)] | units.km / units.s
    

class StellarEvolutionCodeForTesting(object):
    
    def __init__(self, particles=Particles(6)):
        particles.mass = 1 | units.MSun
        particles.radius = range(1, len(particles)+1) | units.RSun
        self.particles = particles
    

class StellarEvolutionCodeWithInternalStructureForTesting(object):
    
    def __init__(self, particles=Particles(6)):
        particles.mass = 1 | units.MSun
        particles.radius = range(1, len(particles)+1) | units.RSun
        particles.type = "native star" | units.string
        self.particles = particles
    
    def new_particle_from_model(self, internal_structure, current_age, key=None):
        tmp_star = Particle(key=key)
        tmp_star.mass = internal_structure["mass"]
        tmp_star.radius = internal_structure["radius"]
        tmp_star.type = "new particle from model" | units.string
        return self.particles.add_particle(tmp_star)
    


class TestCollisionHandler(TestCase):
    
    def test1(self):
        print "Test CollisionHandler with collision code class (creates new instance for each collision)"
        colliders = Particles(2)
        handler = CollisionHandler(CollisionCodeForTesting)
        
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertTrue(isinstance(result, Particles))
        self.assertEqual(result.mass, 1 | units.kg)
        
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertEqual(result.mass, 1 | units.kg)
    
    def test2(self):
        print "Test CollisionHandler with collision code instance (same instance for each collision)"
        colliders = Particles(2)
        handler = CollisionHandler(CollisionCodeForTesting())
        
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertTrue(isinstance(result, Particles))
        self.assertEqual(result.mass, 1 | units.kg)
        
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertEqual(result.mass, 2 | units.kg)
    
    def test3(self):
        print "Test CollisionHandler with collision code class, arguments and parameters"
        colliders = Particles(2)
        handler = CollisionHandler(
            CollisionCodeForTesting, 
            collision_code_arguments=dict(next_mass=5|units.kg),
            collision_code_parameters=dict(mass_unit=units.g)
        )
        
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertTrue(isinstance(result, Particles))
        self.assertEqual(result.mass, 5 | units.kg)
        self.assertTrue(result.mass.unit is units.g)
        
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertEqual(result.mass, 5 | units.kg)
        self.assertTrue(result.mass.unit is units.g)
        
        handler.collision_code_arguments = dict(next_mass=42|units.kg)
        handler.collision_code_parameters = dict(mass_unit=units.MSun)
        result = handler.handle_collision(colliders[0], colliders[1])
        self.assertEqual(result.mass, 42 | units.kg)
        self.assertTrue(result.mass.unit is units.MSun)
    
    def test4(self):
        print "Test handle_collisions"
        colliders = Particles(8)
        handler = CollisionHandler(CollisionCodeForTesting)
        
        result = handler.handle_collisions(colliders[:4], colliders[4:])
        self.assertTrue(isinstance(result, Particles))
        self.assertEqual(len(result), 4)
        self.assertEqual(result.mass, [1, 1, 1, 1] | units.kg)
    
    def test5(self):
        print "Test CollisionHandler with gravity code"
        gravity = GravityCodeForTesting()
        self.assertEqual(len(gravity.particles), 6)
        
        handler = CollisionHandler(CollisionCodeForTesting, gravity_code=gravity)
        merged = handler.handle_collisions(gravity.particles[::2], gravity.particles[1::2])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged.mass, [1, 1, 1] | units.kg)
        
        self.assertEqual(len(gravity.particles), 3)
        self.assertEqual(gravity.particles.mass, [1, 1, 1] | units.kg)
        self.assertEqual(gravity.particles.radius, [3, 3, 3] | units.RSun)
        self.assertAlmostEqual(gravity.particles.position, 
            [[0.5, 1.0, 1.5], [2.5, 5.0, 7.5], [4.5, 9.0, 13.5]] | units.AU)
        self.assertAlmostEqual(gravity.particles.velocity, 
            [[0.5, 0.5, 0.5], [2.5, 6.5, 17.5], [4.5, 20.5, 94.5]] | units.km / units.s)
    
    def test6(self):
        print "Test CollisionHandler with stellar evolution code, type I"
        stellar_evolution = StellarEvolutionCodeForTesting()
        self.assertEqual(len(stellar_evolution.particles), 6)
        
        collision_code = CollisionCodeForTesting()
        collision_code.stellar_evolution_code_required = True
        
        self.assertRaises(AmuseException, CollisionHandler, collision_code, expected_message=
             "CollisionCodeForTesting requires a stellar evolution code: "
             "CollisionHandler(..., stellar_evolution_code=x)")
        
        handler = CollisionHandler(collision_code, stellar_evolution_code=stellar_evolution)
        merged = handler.handle_collisions(stellar_evolution.particles[::2], stellar_evolution.particles[1::2])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged.mass, [1, 2, 3] | units.kg)
        
        self.assertEqual(len(stellar_evolution.particles), 3)
        self.assertEqual(stellar_evolution.particles.mass, [1, 2, 3] | units.kg)
        self.assertEqual(stellar_evolution.particles.radius, [3, 7, 11] | units.RSun)
    
    def test7(self):
        print "Test CollisionHandler with stellar evolution code, type II"
        stellar_evolution = StellarEvolutionCodeWithInternalStructureForTesting()
        self.assertEqual(len(stellar_evolution.particles), 6)
        self.assertEqual(stellar_evolution.particles.type, ["native star"]*6  | units.string)
        
        collision_code = CollisionCodeForTesting()
        collision_code.stellar_evolution_code_required = True
        
        self.assertRaises(AmuseException, CollisionHandler, collision_code, expected_message=
             "CollisionCodeForTesting requires a stellar evolution code: "
             "CollisionHandler(..., stellar_evolution_code=x)")
        
        handler = CollisionHandler(collision_code, stellar_evolution_code=stellar_evolution)
        merged = handler.handle_collisions(stellar_evolution.particles[::2], stellar_evolution.particles[1::2])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged.mass, [1, 2, 3] | units.kg)
        
        self.assertEqual(len(stellar_evolution.particles), 3)
        self.assertEqual(stellar_evolution.particles.mass, [1, 2, 3] | units.kg)
        self.assertEqual(stellar_evolution.particles.radius, [3, 7, 11] | units.RSun)
        self.assertEqual(stellar_evolution.particles.type, ["new particle from model"]*3 | units.string)
    
    def test8(self):
        print "Test CollisionHandler with stellar evolution and gravity code"
        gravity = GravityCodeForTesting()
        self.assertEqual(len(gravity.particles), 6)
        
        stellar_evolution = StellarEvolutionCodeForTesting(particles=gravity.particles.copy())
        self.assertEqual(len(stellar_evolution.particles), 6)
        
        collision_code = CollisionCodeForTesting()
        collision_code.stellar_evolution_code_required = True
        
        self.assertRaises(AmuseException, CollisionHandler, collision_code, expected_message=
             "CollisionCodeForTesting requires a stellar evolution code: "
             "CollisionHandler(..., stellar_evolution_code=x)")
        
        handler = CollisionHandler(
            collision_code, 
            stellar_evolution_code=stellar_evolution, 
            gravity_code=gravity
        )
        merged = handler.handle_collisions(stellar_evolution.particles[::2], stellar_evolution.particles[1::2])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged.mass, [1, 2, 3] | units.kg)
        
        self.assertEqual(len(stellar_evolution.particles), 3)
        self.assertEqual(stellar_evolution.particles.mass, [1, 2, 3] | units.kg)
        self.assertEqual(stellar_evolution.particles.radius, [3, 7, 11] | units.RSun)
        
        self.assertEqual(len(gravity.particles), 3)
        self.assertEqual(gravity.particles.mass, [1, 2, 3] | units.kg)
        self.assertEqual(gravity.particles.radius, [3, 7, 11] | units.RSun)
        self.assertAlmostEqual(gravity.particles.position, 
            [[0.5, 1.0, 1.5], [2.5, 5.0, 7.5], [4.5, 9.0, 13.5]] | units.AU)
        self.assertAlmostEqual(gravity.particles.velocity, 
            [[0.5, 0.5, 0.5], [2.5, 6.5, 17.5], [4.5, 20.5, 94.5]] | units.km / units.s)
    
    def test9(self):
        print "Test CollisionHandler with gravity code and StickySpheres collision code"
        gravity = GravityCodeForTesting()
        self.assertEqual(len(gravity.particles), 6)
        gravity.particles.mass = [1, 1, 2, 2, 3, 3,] | units.MSun
        gravity.particles.radius = range(101, 107) | units.RSun
        
        collision_code = StickySpheres(mass_loss=0.1)
        
        handler = CollisionHandler(
            collision_code, 
            gravity_code=gravity
        )
        merged = handler.handle_collisions(gravity.particles[::2], gravity.particles[1::2])
        self.assertTrue(isinstance(merged, Particles))
        self.assertEqual(len(merged), 3)
        self.assertEqual(merged.mass, 0.9 * ([2, 4, 6] | units.MSun))
        
        self.assertEqual(len(gravity.particles), 3)
        self.assertEqual(gravity.particles.mass, 0.9 * ([2, 4, 6] | units.MSun))
        self.assertEqual(gravity.particles.radius, [102, 104, 106] | units.RSun)
        self.assertAlmostEqual(gravity.particles.position, 
            [[0.5, 1.0, 1.5], [2.5, 5.0, 7.5], [4.5, 9.0, 13.5]] | units.AU)
        self.assertAlmostEqual(gravity.particles.velocity, 
            [[0.5, 0.5, 0.5], [2.5, 6.5, 17.5], [4.5, 20.5, 94.5]] | units.km / units.s)
    
