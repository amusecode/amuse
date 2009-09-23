import unittest

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import parameters

class TestAttributeParameterDefintions(unittest.TestCase):
    def test1(self):
        x = parameters.ModuleAttributeParameterDefinition("test", "test_name", "a test parameter", units.m, 0.1 | units.m)
        self.assertEqual(x.name,'test_name')
        class TestModule(object):
            test = 123
        o = TestModule()
        value = x.get_value(o)
        self.assertEqual(value.number, 123)
        self.assertTrue(value.unit.has_same_base_as(units.m))
        
    def test2(self):
        x = parameters.ModuleAttributeParameterDefinition("test", "test_name", "a test parameter", nbody_system.length, 0.1 | units.m)
        self.assertEqual(x.name,'test_name')
        class TestModule(object):
            convert_nbody = nbody_system.nbody_to_si(1.0 | units.km, 60 | units.s)
            test = 123
        o = TestModule()
        value = x.get_value(o)
        self.assertAlmostEqual(value.in_(units.m).number, 123000, 3)
        self.assertTrue(value.unit.has_same_base_as(units.m))
        
    def test3(self):
        x = parameters.ModuleAttributeParameterDefinition("test", "test_name", "a test parameter", units.m, 0.1 | units.m)
        self.assertEqual(x.name,'test_name')
        class TestModule(object):
            test = 123
        o = TestModule()
        value = x.set_value(o, 5 | units.km)
        self.assertEqual(o.test, 5000)

        
    def test4(self):
        x = parameters.ModuleAttributeParameterDefinition("test", "test_name", "a test parameter", nbody_system.length, 0.1 | units.m)
        class TestModule(object):
            convert_nbody = nbody_system.nbody_to_si(1.0 | units.km, 60 | units.s)
            test = 123
        o = TestModule()
        value = x.set_value(o, 20 | units.km)
        self.assertAlmostEqual(o.test, 20, 10)
        
    def test5(self):
        definition = parameters.ModuleAttributeParameterDefinition("test", "test_name", "a test parameter", units.m, 0.1 | units.m)
        
        class TestModule(object):
            test = 123
        o = TestModule()
        x = parameters.Parameters([definition], o)
        value = x.test_name
        self.assertEqual(value.number, 123)
        x.test_name = 1 | units.km
        self.assertEqual(o.test, 1000)
        

class TestMethodParameterDefintions(unittest.TestCase):
    def test1(self):
        x = parameters.ModuleMethodParameterDefinition("get_test","set_test", "test_name", "a test parameter", units.m, 0.1 | units.m)
        
        class TestModule(object):
            def get_test(self):
                return 123
                
        o = TestModule()
        value = x.get_value(o)
        self.assertEqual(value.number, 123)
        self.assertTrue(value.unit.has_same_base_as(units.m))
        
    def test2(self):
        x = parameters.ModuleMethodParameterDefinition("get_test","set_test", "test_name", "a test parameter", units.m, 0.1 | units.m)
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        x.set_value(o, 10|units.m)
        self.assertEqual(o.x, 10)
        value = x.get_value(o)
        self.assertTrue(value.number, 10)
        
        
    def test3(self):
        x = parameters.ModuleMethodParameterDefinition("get_test","set_test", "test_name", "a test parameter", units.no_unit, 0.1 | units.no_unit)
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        x.set_value(o, 10|units.no_unit)
        self.assertEqual(o.x, 10)
        value = x.get_value(o)
        self.assertTrue(value.number, 10)
        
        
