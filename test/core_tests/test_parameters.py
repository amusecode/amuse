import warnings

from amuse.test import amusetest
from amuse.support.exceptions import AmuseException, AmuseWarning
from amuse.units import nbody_system, generic_unit_system, generic_unit_converter
from amuse.units import units
from amuse.datamodel import parameters
from amuse.support.interface import HandleParameters

class BaseTestModule(object):
    def before_get_parameter(self):
        return
        
    def before_set_parameter(self):
        return
        
class TestMethodParameterDefintions(amusetest.TestCase):
    def test1(self):
        class TestModule(BaseTestModule):
            def get_test(self):
                return 123 | units.m
    
        o = TestModule()
        set = parameters.Parameters([parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | units.m)], o)
        x = set.get_parameter("test_name")
        value = x.get_value()
        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)

    def test2(self):
        definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | units.m)
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        x = set.get_parameter("test_name")
        x.set_value(10|units.m)
        self.assertEqual(o.x, 10|units.m)
        value = x.get_value()
        self.assertEqual(value, 10|units.m)


    def test3(self):
        definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | units.no_unit)
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
    
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        x = set.get_parameter("test_name")
        x.set_value(10|units.none)
        self.assertEqual(o.x, 10|units.none)
        value = x.get_value()
        self.assertEqual(value, 10)

    def test4(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | units.m
        )

        class TestModule(BaseTestModule):
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

        self.assertEquals(1 | units.km, instance.parameters.test_name)
        self.assertEquals(1000 | units.m, instance.x)


    def test5(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            None,
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | units.m
        )

        class TestModule(BaseTestModule):
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

        self.assertEquals(1 | units.km, instance.parameters.test_name)
        self.assertEquals(1000 | units.m, instance.x)



    def test6(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            "bla"
        )

        class TestModule(BaseTestModule):
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


        instance.parameters.test_name = "bla"

        self.assertEquals("bla", instance.x)

        instance.parameters.test_name = "bla"
        self.assertEquals("bla", instance.x )

  
    def test8(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            11.0 | units.m
        )

        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value


        instance = TestModule()

        p = parameters.Parameters([parameter_definition], instance)

        p.set_defaults()

        self.assertEquals(11.0 | units.m, instance.x)


    def test9(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            11.0 | units.m
        )

        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value


        instance = TestModule()

        p = parameters.Parameters([parameter_definition], instance)

        self.assertRaises(AmuseException, lambda: p.unknown, 
            expected_message = "tried to get unknown parameter 'unknown' for a 'TestModule' object")

        with warnings.catch_warnings(record=True) as w:

            p.unknown = 10 | units.m
            self.assertEquals(len(w), 1)
            self.assertEquals("tried to set unknown parameter 'unknown' for a 'TestModule' object", str(w[-1].message))

    def test10(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            11.0 | units.m
        )
    
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
    
    
        instance = TestModule()
    
        p = parameters.Parameters([parameter_definition], instance)
        instance.x = 1 | units.m
        self.assertEquals(p.test_name, 1 | units.m)
        
        def try_set_read_only_parameter(parameter_set):
            parameter_set.test_name = 2 | units.m
        
        self.assertRaises(AmuseException, try_set_read_only_parameter, p, 
            expected_message = "Could not set value for parameter 'test_name' of a 'TestModule' object, parameter is read-only")


    def test11(self):
        parameter_definition1 = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            11.0 | units.m
        )
        
        parameter_definition2 = parameters.ModuleMethodParameterDefinition(
            "get_test1",
            "set_test1",
            "test_name2",
            "a test parameter",
            12.0 | units.m
        )
    
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
            def get_test1(self):
                return self.y
            def set_test1(self, value):
                self.y = value
    
    
        instance = TestModule()
    
        p = parameters.Parameters([parameter_definition1, parameter_definition2], instance)
        instance.x = 1 | units.m
        instance.y = 2 | units.m
        self.assertEquals(p.test_name, 1 | units.m)
        self.assertEquals(p.test_name2, 2 | units.m)
        
        p.test_name = 20 | units.m
        p.send_not_set_parameters_to_code()
        
        self.assertEquals(instance.x, 20 | units.m)
        self.assertEquals(instance.y, 12 | units.m)
        



class TestParameters(amusetest.TestCase):
    def test1(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            11.0 | units.m
        )

        class TestModule(BaseTestModule):
            x = 123 | units.m

            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value


        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)

        value = x.test_name

        self.assertTrue(value.unit.has_same_base_as(units.m))
        self.assertEqual(value.value_in(units.m), 123)

    def test2(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            11.0 | nbody_system.length
        )

        class TestModule(BaseTestModule):
            x = 123 | nbody_system.length

            def get_test(self):
                return self.x
                
            def set_test(self, value):
                self.x = value

        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)


        self.assertEqual(x.test_name, 123 | nbody_system.length)

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)

        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_generic()
            )

        self.assertAlmostEquals(y.test_name.value_in(units.m), 246.0, 6)
        y.test_name = 500 | units.m


        self.assertAlmostEquals(y.test_name.value_in(units.m), 500.0, 6)
        print x.test_name, o.x
        self.assertAlmostEquals(x.test_name.value_in(nbody_system.length), 250.0, 6)
        self.assertAlmostEquals(o.x, 250.0 | nbody_system.length, 6)


    def test3(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            11.0 | nbody_system.length
        )

        class TestModule(BaseTestModule):
            x = 123 | units.m

            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value

        o = TestModule()
        x = parameters.new_parameters_instance_with_docs([parameter_definition], o)

        self.assertTrue("test_name" in x.__doc__)
        self.assertTrue("a test parameter" in x.__doc__)
        self.assertTrue("default" in x.__doc__)
        self.assertTrue("11.0 length" in x.__doc__)

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.new_parameters_with_units_converted_instance_with_docs(
                x,
                convert_nbody.as_converter_from_si_to_generic()
            )

        self.assertTrue("test_name" in y.__doc__)
        self.assertTrue("a test parameter" in y.__doc__)
        self.assertTrue("default" in y.__doc__)
        self.assertTrue("22.0 m" in y.__doc__)

    def test3b(self):
        # Same test as test3, but testing on the class, not instance
        # This makes sure the python 'help' functionality works on parameters
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            11.0 | nbody_system.length
        )

        class TestModule(BaseTestModule):
            x = 123 | units.m

            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value

        o = TestModule()
        x = parameters.new_parameters_instance_with_docs([parameter_definition], o)

        self.assertTrue("test_name" in x.__class__.__doc__)
        self.assertTrue("a test parameter" in x.__class__.__doc__)
        self.assertTrue("default" in x.__class__.__doc__)
        self.assertTrue("11.0 length" in x.__class__.__doc__)

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.new_parameters_with_units_converted_instance_with_docs(
                x,
                convert_nbody.as_converter_from_si_to_generic()
            )

        self.assertTrue("test_name" in y.__class__.__doc__)
        self.assertTrue("a test parameter" in y.__class__.__doc__)
        self.assertTrue("default" in y.__class__.__doc__)
        self.assertTrue("22.0 m" in y.__class__.__doc__)

    def test4(self):
        parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            None,
            "test_name",
            "a test parameter",
            11.0 | nbody_system.length
        )

        class TestModule(BaseTestModule):
            x = 123.0 | nbody_system.length

            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value

        o = TestModule()
        x = parameters.Parameters([parameter_definition], o)


        self.assertTrue("test_name" in str(x))
        self.assertTrue("123.0 length" in str(x))

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_generic()
            )
        self.assertTrue("test_name" in str(y))
        self.assertTrue("246.0 m" in str(y))

    def test5(self):
        print "Test 5: testing mixed nbody and physical units"
        phys_parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "phys_test_name",
            "a test parameter with physical units",
            11.0 | units.m
        )
        nbody_parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_test1",
            "set_test1",
            "nbody_test_name",
            "a test parameter with nbody units",
            11.0 | nbody_system.length
        )

        class TestModule(BaseTestModule):
            x = 123.0 | units.m
            y = 123.0 | nbody_system.length

            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
            def get_test1(self):
                return self.y
            def set_test1(self, value):
                self.y = value

        o = TestModule()
        x = parameters.Parameters([phys_parameter_definition, nbody_parameter_definition], o)

        self.assertTrue("nbody_test_name" in str(x))
        self.assertTrue("123.0 length" in str(x))
        self.assertTrue("phys_test_name" in str(x))
        self.assertTrue("123.0 m" in str(x))

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_generic()
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
        nbody_parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_nbody",
            None,
            "nbody_par_name",
            "a test parameter with nbody units",
            11.0 | nbody_system.length
        )
        string_parameter_definition = parameters.ModuleMethodParameterDefinition(
            "get_string",
            None,
            "string_par_name",
            "a test parameter with string units",
            "test string"
        )

        class TestModule(BaseTestModule):
            x = 123.0 | nbody_system.length

            def get_nbody(self):
                return self.x
            def get_string(self):
                return str(10 * self.x.number )

        o = TestModule()
        x = parameters.Parameters([string_parameter_definition, nbody_parameter_definition], o)

        self.assertTrue("nbody_par_name" in str(x))
        self.assertTrue("123.0 length" in str(x))
        self.assertTrue("string_par_name" in str(x))
        self.assertTrue("1230.0" in str(x))

        convert_nbody = nbody_system.nbody_to_si(2.0 | units.m, 4.0 | units.kg)
        y = parameters.ParametersWithUnitsConverted(
                x,
                convert_nbody.as_converter_from_si_to_generic()
            )
        self.assertEquals(getattr(y,"string_par_name"), "1230.0")
        self.assertAlmostEquals(getattr(y,"nbody_par_name"), 246.0 | units.m)

    def test7(self):
        parameter_definition1 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg1",
            "test_par1",
            "a test parameter (1)",
            11.0 | units.m
        )
    
        parameter_definition2 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg2",
            "test_par2",
            "a test parameter (2)",
            12.0 | units.m
        )
    
        class TestModule(BaseTestModule):
            x = 123 | units.m
            y = 456 | units.m
    
            def initialize_vars(self, arg1, arg2):
                self.x = arg1
                self.y = arg2
    
    
        o = TestModule()
        x = parameters.Parameters([parameter_definition1, parameter_definition2], o)
        x.test_par1 = 20 | units.m
        print x.test_par1
        self.assertEquals(x.test_par1, 20 | units.m)
        self.assertEquals(x.test_par2, 12 | units.m)
        self.assertEquals(o.x, 123 | units.m)
        self.assertEquals(o.y, 456 | units.m)
        x.send_cached_parameters_to_code()
        self.assertEquals(o.x, 20 | units.m)
        self.assertEquals(o.y, 12 | units.m)
    
    

    def test8(self):
        parameter_definition1 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg1",
            "test_par1",
            "a test parameter (1)",
            11.0 | units.m
        )
    
        parameter_definition2 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars",
            "arg2",
            "test_par2",
            "a test parameter (2)",
            12.0 | units.m
        )
    
        parameter_definition3 = parameters.ModuleCachingParameterDefinition(
            "initialize_vars2",
            "arg1",
            "test_par3",
            "a test parameter (3)",
            14.0 | units.m
        )
    
        class TestModule(BaseTestModule):
            x = 123 | units.m
            y = 456 | units.m
            z = 100 | units.m
    
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
        self.assertEquals(o.x, 11 | units.m)
        self.assertEquals(o.y, 12 | units.m)
        self.assertEquals(o.z, 14 | units.m)
    
    

    def test9(self):
        parameter_definition1 = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            11.0 | units.m
        )
        
        parameter_definition2 = parameters.ModuleMethodParameterDefinition(
            "get_test1",
            "set_test1",
            "test_name2",
            "a test parameter",
            12.0 | units.m
        )
    
        paramer_definition3 = parameters.VectorParameterDefinition(
            "test_vector",
            "vector of parameters",
            ["test_name", "test_name2"],
            [11.0, 12.0] | units.m
        )
        
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
            def get_test1(self):
                return self.y
            def set_test1(self, value):
                self.y = value
    
    
        instance = TestModule()
        instance.x = 1 | units.m
        instance.y = 2 | units.m
        
        p = parameters.Parameters([parameter_definition1, parameter_definition2, paramer_definition3], instance)
       
        self.assertEquals(p.test_vector, (1,2) | units.m)
        p.test_vector = (3,4) | units.m
        self.assertEquals(instance.x, 3 | units.m)
        self.assertEquals(instance.y, 4 | units.m)
    
    def test10(self):
        print "Testing ParametersWithUnitsConverted on vector parameters"
        definitions = []
        for par_name in ["length_x", "length_y", "length_z"]:
            definitions.append(parameters.ModuleMethodParameterDefinition(
                "get_"+par_name,
                "set_"+par_name,
                par_name,
                "a test parameter",
                10.0 | generic_unit_system.length
            ))
        
        definitions.append(parameters.VectorParameterDefinition(
            "mesh_length",
            "length of the model in the x, y and z directions",
            ("length_x", "length_y", "length_z"),
            [10, 10, 10] | generic_unit_system.length
        ))
        
        class TestModule(BaseTestModule):
            x = 123.0 | generic_unit_system.length
            y = 456.0 | generic_unit_system.length
            z = 789.0 | generic_unit_system.length
            
            def get_length_x(self):
                return self.x
            def set_length_x(self, value):
                self.x = value
            def get_length_y(self):
                return self.y
            def set_length_y(self, value):
                self.y = value
            def get_length_z(self):
                return self.z
            def set_length_z(self, value):
                self.z = value
        
        o = TestModule()
        x = parameters.Parameters(definitions, o)
        
        self.assertTrue("mesh_length" in str(x))
        self.assertTrue("[123.0, 456.0, 789.0] length" in str(x))
        
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(2.0 | units.m, 4.0 | units.kg, 6.0 | units.s)
        y = parameters.ParametersWithUnitsConverted(
                x,
                converter.as_converter_from_si_to_generic()
            )
        self.assertTrue("mesh_length" in str(y))
        self.assertTrue("[246.0, 912.0, 1578.0] m" in str(y))
    
    
    def test11(self):
        print "Testing ParametersWithUnitsConverted on vector parameters, using add_vector_parameter"
        
        class TestModule(BaseTestModule):
            x = 123.0 | generic_unit_system.length
            y = 456.0 | generic_unit_system.length
            z = 789.0 | generic_unit_system.length
            
            def get_length_x(self):
                return self.x
            def set_length_x(self, value):
                self.x = value
            def get_length_y(self):
                return self.y
            def set_length_y(self, value):
                self.y = value
            def get_length_z(self):
                return self.z
            def set_length_z(self, value):
                self.z = value
        
        o = TestModule()
        parameters_handler = HandleParameters(o)
        parameters_handler.add_vector_parameter(
            "mesh_length",
            "length of the model in the x, y and z directions",
            ("length_x", "length_y", "length_z")
        )
        for par_name in ["length_x", "length_y", "length_z"]:
            parameters_handler.add_method_parameter(
                "get_"+par_name, 
                "set_"+par_name,
                par_name, 
                "a test parameter", 
                default_value = 10.0 | generic_unit_system.length,
            )
        
        x = parameters_handler.get_attribute(None, None)
        self.assertTrue("mesh_length" in str(x))
        self.assertTrue("[123.0, 456.0, 789.0] length" in str(x))
        
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(2.0 | units.m, 4.0 | units.kg, 6.0 | units.s)
        y = parameters.ParametersWithUnitsConverted(
                x,
                converter.as_converter_from_si_to_generic()
            )
        self.assertTrue("mesh_length" in str(y))
        self.assertTrue("[246.0, 912.0, 1578.0] m" in str(y))
    
    def test12(self):
        definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            "set_test",
            "test_name",
            "a test parameter",
            0.1 | units.m
        )
        class TestModule(BaseTestModule):
            def get_test(self):
                return self.x
            def set_test(self, value):
                self.x = value
                
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        set.test_name = 10|units.m
        
        self.assertEqual(o.x, 10|units.m)
        self.assertEqual(set.test_name, 10|units.m)
        
        memento = set.copy()
        self.assertEqual(memento.test_name, 10|units.m)
        set.test_name = 20|units.m
        
        self.assertEqual(o.x, 20|units.m)
        self.assertEqual(set.test_name, 20|units.m)
        self.assertEqual(memento.test_name, 10|units.m)
        
        set.reset_from_memento(memento)
        
        self.assertEqual(o.x, 10|units.m)
        self.assertEqual(set.test_name, 10|units.m)
        self.assertEqual(memento.test_name, 10|units.m)
    
    def test13(self):
        definition = parameters.ModuleMethodParameterDefinition(
            "get_test",
            None,
            "test_name",
            "a read-only test parameter",
            0.1 | units.m
        )
        class TestModule(BaseTestModule):
            x = 0.1 | units.m
            def get_test(self):
                return self.x
        
        o = TestModule()
        set = parameters.Parameters([definition,], o)
        
        self.assertRaises(AmuseException, setattr, set, "test_name", 1.0 | units.m,
            expected_message = "Could not set value for parameter 'test_name' of a 'TestModule' object, parameter is read-only")
        
        self.assertEqual(o.x, 0.1|units.m)
        self.assertEqual(set.test_name, 0.1|units.m)
        
        memento = set.copy()
        self.assertEqual(memento.test_name, 0.1|units.m)
        
        set.reset_from_memento(memento)
        self.assertEqual(o.x, 0.1|units.m)
        self.assertEqual(set.test_name, 0.1|units.m)
        
        memento.test_name = 2.0 | units.m
        self.assertEqual(memento.test_name, 2.0|units.m)
        
        with warnings.catch_warnings(record=True) as w:
            set.reset_from_memento(memento)
            self.assertEquals(len(w), 1)
            self.assertEquals("tried to change read-only parameter 'test_name' for a 'TestModule' object", str(w[-1].message))
        
        self.assertEqual(o.x, 0.1|units.m)
        self.assertEqual(set.test_name, 0.1|units.m)
        self.assertEqual(memento.test_name, 2.0|units.m)


