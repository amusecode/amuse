from amuse.test import amusetest
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import parameters

import warnings
from amuse.support import exceptions

class TestAttributeParameterDefintions(amusetest.TestCase):
    def test1(self):
        class TestModule(object):
            def __init__(self):
                self.legacy_interface = self
                self.test = 123
        o = TestModule()
        set = parameters.Parameters([parameters.ModuleAttributeParameterDefinition("test", "test_name", "a test parameter", units.m, 0.1 | units.m)], o)
        x = set.get_parameter("test_name")
        value = x.get_value()
        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)

    def test2(self):
        class TestModule(object):
            def __init__(self):
                self.legacy_interface = self
                self.test = 123
    
        o = TestModule()
        set = parameters.Parameters([parameters.ModuleAttributeParameterDefinition(
            "test",
            "test_name",
            "a test parameter",
            nbody_system.length,
            0.1 | nbody_system.length
        )], o)
        x = set.get_parameter("test_name")
        value = x.get_value()
        self.assertEquals(value.value_in(nbody_system.length), 123.0)
        self.assertTrue(value.unit.has_same_base_as(nbody_system.length))

    def test3(self):
        class TestModule(object):
            def __init__(self):
                self.legacy_interface = self
                self.test = 123
        o = TestModule()
        set = parameters.Parameters([parameters.ModuleAttributeParameterDefinition(
            "test",
            "test_name",
            "a test parameter",
            units.m, 0.1 | units.m
        )], o)
        x = set.get_parameter("test_name")
        value = x.set_value(5 | units.km)
        self.assertEqual(o.test, 5000)


    def test4(self):
        class TestModule(object):
            def __init__(self):
                self.legacy_interface = self
                self.test = 123
    
        o = TestModule()
        set = parameters.Parameters([parameters.ModuleAttributeParameterDefinition(
            "test",
            "test_name",
            "a test parameter",
            nbody_system.length,
            0.1 | nbody_system.length
        )], o)
        x = set.get_parameter("test_name")
        value = x.set_value(2 | (10.0 * nbody_system.length) )
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
            def __init__(self):
                self.legacy_interface = self
                self.test = 123
    
        o = TestModule()
        x = parameters.Parameters([definition], o)
        value = x.test_name
        self.assertEqual(value.value_in(units.m), 123)
        x.test_name = 1 | units.km
        self.assertEqual(o.test, 1000)


class TestMethodParameterDefintions(amusetest.TestCase):
    def test1(self):
        class TestModule(object):
            def get_test(self):
                return 123,0
    
        o = TestModule()
        set = parameters.Parameters([parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            0.1 | units.m)], o)
        x = set.get_parameter("test_name")
        value = x.get_value()
        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)

    def test2(self):
        definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            0.1 | units.m)
        class TestModule(object):
            def get_test(self):
                return self.x,0
            def set_test(self, value):
                self.x = value
                return 0
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        x = set.get_parameter("test_name")
        x.set_value(10|units.m)
        self.assertEqual(o.x, 10)
        value = x.get_value()
        self.assertTrue(value.value_in(units.m), 10)


    def test3(self):
        definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.no_unit,
            0.1 | units.no_unit)
        class TestModule(object):
            def get_test(self):
                return self.x, 0
            def set_test(self, value):
                self.x = value
                return 0
    
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        x = set.get_parameter("test_name")
        x.set_value(10|units.none)
        self.assertEqual(o.x, 10)
        value = x.get_value()
        self.assertTrue(value.value_in(units.none), 10)

    def test4(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            0.1 | units.m
        )

        class TestModule(object):
            def get_test(self):
                return self.x, 0
            def set_test(self, value):
                self.x = value
                return 0

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
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            None,
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            0.1 | units.m
        )

        class TestModule(object):
            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0

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
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.string,
            "bla" | units.string
        )

        class TestModule(object):
            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0

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
        definition = parameters.ModuleMethodParameterDefinition_Next(
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
        set = parameters.Parameters([definition,], o)
        x = set.get_parameter("test_name")
        o.errorcode = 0
        x.set_value(10 | units.m)
        self.assertEqual(o.x, 10)
        value = x.get_value()
        self.assertTrue(value.value_in(units.m), 10)
    
        o.errorcode = -1
        try:
            x.set_value(10 | units.m)
            self.fail("Setting the value should result in an exception as the errorcode is set")
        except parameters.ParameterException as ex:
            self.assertEquals(-1, ex.errorcode)
            self.assertEquals("test_name", ex.parameter_name)
            self.assertEquals("Could not set value for parameter 'test_name' of a 'TestModule' object, got errorcode <-1>", str(ex))
    
        o.errorcode = -2
        try:
            x.get_value()
            self.fail("Gettting the value should result in an exception as the errorcode is set")
        except parameters.ParameterException as ex:
            self.assertEquals(-2, ex.errorcode)
            self.assertEquals("test_name", ex.parameter_name)
            self.assertEquals("Could not get value for parameter 'test_name' of a 'TestModule' object, got errorcode <-2>", str(ex))

    def test8(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            11.0 | units.m
        )

        class TestModule(object):
            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0


        instance = TestModule()

        p = parameters.Parameters([parameter_definition], instance)

        p.set_defaults()

        self.assertEquals(11.0 , instance.x)


    def test9(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            11.0 | units.m
        )

        class TestModule(object):
            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0


        instance = TestModule()

        p = parameters.Parameters([parameter_definition], instance)

        try:
            p.unknown
            self.fail("Gettting the value of an unknown parameter should result in an exception")
        except Exception as ex:
            self.assertEquals("tried to get unknown parameter 'unknown' for a 'TestModule' object", str(ex))

        with warnings.catch_warnings(record=True) as w:

            p.unknown = 10 | units.m
            self.assertEquals(len(w), 1)
            self.assertEquals("tried to set unknown parameter 'unknown' for a 'TestModule' object", str(w[-1].message))

    def test10(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            units.m,
            11.0 | units.m
        )
    
        class TestModule(object):
            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0
    
    
        instance = TestModule()
    
        p = parameters.Parameters([parameter_definition], instance)
        instance.x = 1
        self.assertEquals(p.test_name , 1 | units.m)
    
        try:
            p.test_name = 2 | units.m
        except exceptions.CoreException, e:
            self.assertEquals("Could not set value for parameter 'test_name' of a 'TestModule' object, parameter is read-only", str(e))
        else:
            self.fail("Should raise readonly exception")

class TestParameters(amusetest.TestCase):
    def test1(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            units.m,
            11.0 | units.m
        )

        class TestModule(object):
            x = 123

            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0


        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)

        value = x.test_name

        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)

    def test2(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            nbody_system.length,
            11.0 | nbody_system.length
        )

        class TestModule(object):
            x = 123

            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0

        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)


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
        self.assertAlmostEquals(o.x, 250.0, 6)


    def test3(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            nbody_system.length,
            11.0 | nbody_system.length
        )

        class TestModule(object):
            x = 123

            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0

        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)


        self.assertTrue("test_name" in x.__doc__)
        self.assertTrue("a test parameter" in x.__doc__)
        self.assertTrue("default" in x.__doc__)
        self.assertTrue("11.0 nbody length" in x.__doc__)

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_nbody()
            )

        self.assertTrue("test_name" in y.__doc__)
        self.assertTrue("a test parameter" in y.__doc__)
        self.assertTrue("default" in y.__doc__)
        self.assertTrue("22.0 m" in y.__doc__)

    def test4(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            nbody_system.length,
            11.0 | nbody_system.length
        )

        class TestModule(object):
            x = 123.0

            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0

        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)


        self.assertTrue("test_name" in str(x))
        self.assertTrue("123.0 nbody length" in str(x))

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_nbody()
            )
        self.assertTrue("test_name" in str(y))
        self.assertTrue("246.0 m" in str(y))

    def test5(self):
        print "Test 5: testing mixed nbody and physical units"
        phys_parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "phys_test_name",
            "a test parameter with physical units",
            units.m,
            11.0 | units.m
        )
        nbody_parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_test",
            "set_test",
            "nbody_test_name",
            "a test parameter with nbody units",
            nbody_system.length,
            11.0 | nbody_system.length
        )

        class TestModule(object):
            x = 123.0

            def get_test(self):
                return (self.x,0)
            def set_test(self, value):
                self.x = value
                return 0

        o = TestModule()
        x = parameters.Parameters([phys_parameter_definition, nbody_parameter_definition], o)

        self.assertTrue("nbody_test_name" in str(x))
        self.assertTrue("123.0 nbody length" in str(x))
        self.assertTrue("phys_test_name" in str(x))
        self.assertTrue("123.0 m" in str(x))

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_nbody()
            )
        self.assertEquals(getattr(y,"phys_test_name"), 123.0 | units.m)
        self.assertAlmostEquals(getattr(y,"nbody_test_name"), 246.0 | units.m)
        y.phys_test_name = 1234.0 | units.m
        self.assertEquals(y.phys_test_name, 1234.0 | units.m)
        y.nbody_test_name = 12345.0 | nbody_system.length
        self.assertAlmostEquals(y.nbody_test_name, 24690.0 | units.m)
        y.nbody_test_name = 12345.0 | units.m
        self.assertEquals(y.nbody_test_name, 12345.0 | units.m)

    def test6(self):
        print "Test 5: testing mixed nbody and string units"
        nbody_parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_nbody",
            None,
            "nbody_par_name",
            "a test parameter with nbody units",
            nbody_system.length,
            11.0 | nbody_system.length
        )
        string_parameter_definition = parameters.ModuleMethodParameterDefinition_Next(
            "get_string",
            None,
            "string_par_name",
            "a test parameter with string units",
            units.string,
            "test string" | units.string
        )

        class TestModule(object):
            x = 123.0

            def get_nbody(self):
                return (self.x,0)
            def get_string(self):
                return (str(10 * self.x),0)

        o = TestModule()
        x = parameters.Parameters([string_parameter_definition, nbody_parameter_definition], o)


        self.assertTrue("nbody_par_name" in str(x))
        self.assertTrue("123.0 nbody length" in str(x))
        self.assertTrue("string_par_name" in str(x))
        self.assertTrue("1230.0" in str(x))

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_nbody()
            )
        self.assertEquals(getattr(y,"string_par_name"), "1230.0" | units.string)
        self.assertAlmostEquals(getattr(y,"nbody_par_name"), 246.0 | units.m)

    def test7(self):
        parameter_definition1 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg1",
            "test_par1",
            "a test parameter (1)",
            units.m,
            11.0 | units.m
        )
    
        parameter_definition2 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg2",
            "test_par2",
            "a test parameter (2)",
            units.m,
            12.0 | units.m
        )
    
        class TestModule(object):
            x = 123
            y = 456
    
            def initialize_vars(self, arg1, arg2):
                self.x = arg1
                self.y = arg2
                return 0
    
    
        o = TestModule()
        x = parameters.Parameters([parameter_definition1, parameter_definition2], o)
        x.test_par1 = 20 | units.m
        print x.test_par2
        self.assertEquals(o.x, 123)
        self.assertEquals(o.y, 456)
        x.send_cached_parameters_to_code()
        self.assertEquals(o.x, 20)
        self.assertEquals(o.y, 12)
    
    

    def test8(self):
        parameter_definition1 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg1",
            "test_par1",
            "a test parameter (1)",
            units.m,
            11.0 | units.m
        )
    
        parameter_definition2 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg2",
            "test_par2",
            "a test parameter (2)",
            units.m,
            12.0 | units.m
        )
    
        parameter_definition3 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars2",
            "arg1",
            "test_par3",
            "a test parameter (3)",
            units.m,
            14.0 | units.m
        )
    
        class TestModule(object):
            x = 123
            y = 456
            z = 100
    
            def initialize_vars(self, arg1, arg2):
                self.x = arg1
                self.y = arg2
                return 0
    
            def initialize_vars2(self, arg1):
                self.z = arg1
                return 0
    
    
        o = TestModule()
        x = parameters.Parameters([parameter_definition1, parameter_definition2, parameter_definition3], o)
        
        x.send_cached_parameters_to_code()
        self.assertEquals(o.x, 11)
        self.assertEquals(o.y, 12)
        self.assertEquals(o.z, 14)
    
    
