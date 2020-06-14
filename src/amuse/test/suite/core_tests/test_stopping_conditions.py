from amuse.test import amusetest

from amuse.support.exceptions import AmuseException
from amuse.community.interface.stopping_conditions import StoppingConditions
from amuse import datamodel
from amuse.units import units
from amuse.support import interface


class TestStoppingCondition(amusetest.TestCase):
    

    def test1(self):
        
        class AllEnabled(object):
            
            def is_stopping_condition_enabled(self, sc_type):
                return 1
            
            def has_stopping_condition(self, sc_type):
                return 1
        
        instance = StoppingConditions(AllEnabled())
        self.assertTrue(instance.collision_detection.is_supported())
        self.assertTrue(instance.collision_detection.is_enabled())
        self.assertTrue(instance.escaper_detection.is_supported())
        self.assertTrue(instance.escaper_detection.is_enabled())
        self.assertTrue(instance.timeout_detection.is_supported())
        self.assertTrue(instance.timeout_detection.is_enabled())
    
    def test2(self):
        
        class OneEnabled(object):
            
            def is_stopping_condition_enabled(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def has_stopping_condition(self, sc_type):
                return 1 if sc_type == 0 else 0
        
        instance = StoppingConditions(OneEnabled())
        self.assertTrue(instance.collision_detection.is_supported())
        self.assertTrue(instance.collision_detection.is_enabled())
        self.assertFalse(instance.escaper_detection.is_supported())
        self.assertFalse(instance.escaper_detection.is_enabled())
        self.assertFalse(instance.timeout_detection.is_supported())
        self.assertFalse(instance.timeout_detection.is_enabled())
        
    
    def test3(self):
        
        class OneSettable(object):
            is_enabled = 0
            
            def is_stopping_condition_enabled(self, sc_type):
                return self.is_enabled if sc_type == 0 else 0
            
            def has_stopping_condition(self, sc_type):
                return 1 if sc_type == 0 else 0
        
            def enable_stopping_condition(self, sc_type):
                if sc_type == 0:
                    self.is_enabled = 1
        
        instance = StoppingConditions(OneSettable())
        self.assertTrue(instance.collision_detection.is_supported())
        self.assertFalse(instance.collision_detection.is_enabled())
        instance.collision_detection.enable()
        self.assertTrue(instance.collision_detection.is_enabled())
    
    
    
    def test4(self):
        
        class OneSettable(object):
            is_enabled = 0
            
            def is_stopping_condition_enabled(self, sc_type):
                return self.is_enabled if sc_type == 0 else 0
            
            def has_stopping_condition(self, sc_type):
                return 1 if sc_type == 0 else 0
        
            def enable_stopping_condition(self, sc_type):
                if sc_type == 0:
                    self.is_enabled = 1
                    
            def disable_stopping_condition(self, sc_type):
                if sc_type == 0:
                    self.is_enabled = 0
        
        instance = StoppingConditions(OneSettable())
        self.assertTrue(instance.collision_detection.is_supported())
        self.assertFalse(instance.collision_detection.is_enabled())
        instance.collision_detection.enable()
        self.assertTrue(instance.collision_detection.is_enabled())
        instance.collision_detection.disable()
        self.assertFalse(instance.collision_detection.is_enabled())
        
    def test5(self):
        
        class OneEnabled(object):
            
            def is_stopping_condition_enabled(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def has_stopping_condition(self, sc_type):
                return 1 if sc_type == 0 else 0
        
        instance = StoppingConditions(OneEnabled())
        self.assertFalse(instance.escaper_detection.is_supported())
        self.assertFalse(instance.escaper_detection.is_enabled())
        self.assertRaises(AmuseException, instance.escaper_detection.enable)
        self.assertRaises(AmuseException, instance.escaper_detection.disable)
        
        
    def test6(self):
        
        class Collision(object):
   
            def is_stopping_condition_enabled(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def has_stopping_condition(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def is_stopping_condition_set(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def get_number_of_stopping_conditions_set(self):
                return 1
        
            def get_stopping_condition_info(self, indices):
                return [0],[1]
        
        instance = StoppingConditions(Collision())
        instance.code.particles = datamodel.Particles(3)
        instance.code.particles.mass = (1,2,3) | units.kg
        instance.code.particles.add_function_attribute(
            "get_stopping_condition_particle_index",
            lambda particles, indices, sc_type : particles[indices]
        )
        self.assertTrue(instance.collision_detection.is_set())
        particles = instance.collision_detection.particles(0)
        self.assertEqual(len(particles),1)
        self.assertAlmostRelativeEqual(particles[0].mass, 1|units.kg)
        
    
    def test7(self):
        
        class Collision(object):
   
            def is_stopping_condition_enabled(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def has_stopping_condition(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def is_stopping_condition_set(self, sc_type):
                return 1 if sc_type == 0 else 0
            
            def get_number_of_stopping_conditions_set(self):
                return 1
        
            def get_stopping_condition_info(self, indices):
                return [0],[3]
        
        instance = StoppingConditions(Collision())
        instance.code.particles = datamodel.Particles(3)
        instance.code.particles.mass = (1,2,3) | units.kg
        instance.code.particles.add_function_attribute(
            "get_stopping_condition_particle_index",
            lambda particles, indices, sc_type : particles[indices]
        )
        self.assertTrue(instance.collision_detection.is_set())
        particles = instance.collision_detection.particles(0)
        self.assertEqual(len(particles),1)
        particles = instance.collision_detection.particles(1)
        self.assertEqual(len(particles),1)
        particles = instance.collision_detection.particles(2)
        self.assertEqual(len(particles),1)
        particles = instance.collision_detection.particles(3)
        self.assertEqual(len(particles),0)
