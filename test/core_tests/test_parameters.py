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
        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)
        
    def test2(self):
        x = parameters.ModuleAttributeParameterDefinition(
            "test", 
            "test_name", 
            "a test parameter", 
            nbody_system.length, 
            0.1 | nbody_system.length
        )
        self.assertEqual(x.name,'test_name')
        class TestModule(object):
            test = 123.0
            
        o = TestModule()
        value = x.get_value(o)
        self.assertEquals(value.value_in(nbody_system.length), 123.0)
        self.assertTrue(value.unit.has_same_base_as(nbody_system.length))
        
    def test3(self):
        x = parameters.ModuleAttributeParameterDefinition(
            "test", 
            "test_name", 
            "a test parameter", 
            units.m, 0.1 | units.m
        )
        self.assertEqual(x.name,'test_name')
        class TestModule(object):
            test = 123
        o = TestModule()
        value = x.set_value(o, 5 | units.km)
        self.assertEqual(o.test, 5000)

        
    def test4(self):
        x = parameters.ModuleAttributeParameterDefinition(
            "test", 
            "test_name", 
            "a test parameter", 
            nbody_system.length, 
            0.1 | nbody_system.length
        )
        class TestModule(object):
            test = 123.0
        o = TestModule()
        value = x.set_value(o, 2 | (10.0 * nbody_system.length) )
        self.assertAlmostEqual(o.test, 20.0, 10)
        
    def test5(self):
        definition = parameters.ModuleAttributeParameterDefinition(
            "test", 
            "test_name", 
            "a test parameter", 
            units.m, 
            0.1 | units.m
        )
        
        class TestModule(object):
            test = 123
        o = TestModule()
        x = parameters.Parameters([definition], o)
        value = x.test_name
        self.assertEqual(value.value_in(units.m), 123)
        x.test_name = 1 | units.km
        self.assertEqual(o.test, 1000)
        

class TestMethodParameterDefintions(unittest.TestCase):
    def test1(self):
        x = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test", 
            "test_name", 
            "a test parameter", 
            units.m, 
            0.1 | units.m)
        
        class TestModule(object):
            def get_test(self):
                return 123
                
        o = TestModule()
        value = x.get_value(o)
        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)
        
    def test2(self):
        x = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test", 
            "test_name", 
            "a test parameter", 
            units.m, 
            0.1 | units.m)
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        x.set_value(o, 10|units.m)
        self.assertEqual(o.x, 10)
        value = x.get_value(o)
        self.assertTrue(value.value_in(units.m), 10)
        
        
    def test3(self):
        x = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test", 
            "test_name", 
            "a test parameter", 
            units.no_unit, 
            0.1 | units.no_unit)
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        x.set_value(o, 10|units.none)
        self.assertEqual(o.x, 10)
        value = x.get_value(o)
        self.assertTrue(value.value_in(units.none), 10)
        
    def test4(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter", 
            units.m, 
            0.1 | units.m
        )
        
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
        
        class TestModuleBinding(object):
            parameter_definitions = [parameter_definition]
            
            def __init__(self):
                self.parameters = parameters.Parameters(self.parameter_definitions, self)
        
        class TestInterface(TestModule, TestModuleBinding):
            
            def __init__(self):
                TestModuleBinding.__init__(self)
                
        instance = TestInterface()
        
        self.assertTrue('test_name' in list(instance.parameters.names()))
        
        instance.parameters.test_name = 1 | units.km
        
        self.assertEquals(1 | units.km , instance.parameters.test_name)
        self.assertEquals(1000 , instance.x)
        
    
    def test5(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            None,
            "set_test",
            "test_name",
            "a test parameter", 
            units.m, 
            0.1 | units.m
        )
        
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
        
        class TestModuleBinding(object):
            parameter_definitions = [parameter_definition]
            
            def __init__(self):
                self.parameters = parameters.Parameters(self.parameter_definitions, self)
        
        class TestInterface(TestModule, TestModuleBinding):
            
            def __init__(self):
                TestModuleBinding.__init__(self)
                
        instance = TestInterface()
        
        self.assertTrue('test_name' in list(instance.parameters.names()))
        
        instance.parameters.test_name = 1 | units.km
        
        self.assertEquals(1 | units.km , instance.parameters.test_name)
        self.assertEquals(1000 , instance.x)
        
    
    
    def test6(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter", 
            units.string, 
            "bla" | units.string
        )
        
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
        
        class TestModuleBinding(object):
            parameter_definitions = [parameter_definition]
            
            def __init__(self):
                self.parameters = parameters.Parameters(self.parameter_definitions, self)
        
        class TestInterface(TestModule, TestModuleBinding):
            
            def __init__(self):
                TestModuleBinding.__init__(self)
                
        instance = TestInterface()
        
        
        instance.parameters.test_name = "bla" | units.string
        
        self.assertEquals("bla" , instance.x)
        
        instance.parameters.test_name = "bla" 
        self.assertEquals("bla" , instance.x)
        
    def test7(self):
        x = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test", 
            "test_name", 
            "a test parameter", 
            units.m, 
            0.1 | units.m)
    
        class TestModule(object):
            def get_test(self):
                return (self.x, self.errorcode)
            def set_test(self, value):
                self.x = value
                return self.errorcode
                
        o = TestModule()
        o.errorcode = 0
        x.set_value(o, 10 | units.m)
        self.assertEqual(o.x, 10)
        value = x.get_value(o)
        self.assertTrue(value.value_in(units.m), 10)
        
        o.errorcode = -1
        try:
            x.set_value(o, 10 | units.m)
            self.fail("Setting the value should result in an exception as the errorcode is set")
        except parameters.ParameterException as ex:
            self.assertEquals(-1, ex.errorcode)
            self.assertEquals("test_name", ex.parameter_name)
            self.assertEquals("Could not set value for parameter 'test_name' of a 'TestModule' object, got errorcode <-1>", str(ex))
        
        o.errorcode = -2
        try:
            x.get_value(o)
            self.fail("Gettting the value should result in an exception as the errorcode is set")
        except parameters.ParameterException as ex:
            self.assertEquals(-2, ex.errorcode)
            self.assertEquals("test_name", ex.parameter_name)
            self.assertEquals("Could not get value for parameter 'test_name' of a 'TestModule' object, got errorcode <-2>", str(ex))
        
    def test8(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter", 
            units.m, 
            11.0 | units.m
        )
        
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
        
                
        instance = TestModule()
        
        p = parameters.Parameters([parameter_definition], instance)
        
        p.set_defaults()
        
        self.assertEquals(11.0 , instance.x)
        
    
    def test9(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter", 
            units.m, 
            11.0 | units.m
        )
        
        class TestModule(object):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
        
                
        instance = TestModule()
        
        p = parameters.Parameters([parameter_definition], instance)
        
        try:
            p.unknown
            self.fail("Gettting the value of an unknown parameter should result in an exception")
        except Exception as ex:
            self.assertEquals("tried to get unknown parameter 'unknown' for a 'TestModule' object", str(ex))
            
        try:
            p.unknown = 10 | units.m
            self.fail("Settting the value of an unknown parameter should result in an exception")
        except Exception as ex:
            self.assertEquals("tried to set unknown parameter 'unknown' for a 'TestModule' object", str(ex))
        
        
class TestParameters(unittest.TestCase):
    def test1(self):
        definition = parameters.ModuleAttributeParameterDefinition(
            "test", 
            "test_name", 
            "a test parameter", 
            units.m, 
            0.1 | units.m
        )
        
        
        class TestModule(object):
            test = 123
            
        o = TestModule()
        x = parameters.Parameters([definition], o)
        
        value = x.test_name
        
        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)
   
    def test2(self):
        definition = parameters.ModuleAttributeParameterDefinition(
            "test", 
            "test_name", 
            "a test parameter", 
            nbody_system.length, 
            0.1 | nbody_system.length
        )
        
        
        class TestModule(object):
            test = 123
            
        o = TestModule()
        x = parameters.Parameters([definition], o)
        
        
        self.assertEqual(x.test_name.value_in(nbody_system.length), 123)
        
        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_nbody()
            )
        
        self.assertAlmostEquals(y.test_name.value_in(units.m), 246.0, 6)
        y.test_name = 500 | units.m
        
        
        self.assertAlmostEquals(y.test_name.value_in(units.m), 500.0, 6)
        self.assertAlmostEquals(x.test_name.value_in(nbody_system.length), 250.0, 6)
        self.assertAlmostEquals(o.test, 250.0, 6)
        
        
