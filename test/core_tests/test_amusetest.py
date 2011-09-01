from amuse.test.amusetest import TestCase
from amuse.support.exceptions import AmuseException
from amuse.support.units import units


class TestAmusetestTestCase(TestCase):
    
    def test1(self):
        self.assertEqual(1 | units.m, 1 | units.m)
        self.assertEqual(1 | units.m, 1.0 | units.m)
        self.assertEqual([1,2,3] | units.m, [1,2,3] | units.m)
        self.assertEqual([1,1,1] | units.m, 1 | units.m)
        self.assertEqual(1 | units.m, [1,1,1] | units.m)
        self.assertRaises(self.failureException, self.assertEqual, 1 | units.m, 0 | units.m)
        self.assertRaises(self.failureException, self.assertEqual, 
            1 | units.m, 0 | units.m, expected_message = "1 m != 0 m")
    
    def test2(self):
        self.assertAlmostEqual(1 | units.m, 1 | units.m)
        self.assertAlmostEqual(1 | units.m, 1.000000049 | units.m)
        self.assertAlmostEqual([1,2,3] | units.m, [1,2,3] | units.m)
        self.assertAlmostEqual([1,1,1] | units.m, 1.000000049 | units.m)
        self.assertAlmostEqual(1 | units.m, [1,1,1.000000049] | units.m)
        self.assertAlmostEqual(1 | units.m, 1.00049 | units.m, places = 3)
        self.assertRaises(self.failureException, self.assertAlmostEqual, 1.0 | units.m, 
            1.000000051 | units.m, expected_message = "1.0 m != 1.000000051 m within 7 places")
        self.assertRaises(self.failureException, self.assertAlmostEqual, 
            [1.1, 1.0, 0.9] | units.m, 1.0 | units.m, 
            expected_message = "@0, 1.1 m != 1.0 m within 7 places\n@2, 0.9 m != 1.0 m within 7 places")
    
    def test3(self):
        self.assertAlmostRelativeEqual(1 | units.m, 1 | units.m)
        self.assertAlmostRelativeEqual(1 | units.m, 1.0 + 4.9e-16 | units.m)
        self.assertAlmostRelativeEqual([1,2,3] | units.m, [1,2,3] | units.m)
        self.assertAlmostRelativeEqual([1,1,1] | units.m, 1 | units.m)
        self.assertAlmostRelativeEqual(1 | units.m, [1,1,1] | units.m)
        self.assertAlmostRelativeEqual(100.0 | units.m, 100.1 | units.m, places = 3)
        self.assertRaises(self.failureException, self.assertAlmostRelativeEqual, 
            1.0 | units.m, [0.0, 1.0001, 1.0011] | units.m, places = 3, expected_message = 
            "@0, quantity<1.0 m> != quantity<0.0 m> within 3 places\n"
            "@2, quantity<1.0 m> != quantity<1.0011 m> within 3 places")
    
    def test4(self):
        self.assertIsOfOrder(1 | units.m, 1 | units.m)
        self.assertIsOfOrder(1 | units.m, 3.0 | units.m)
        self.assertIsOfOrder([1,2,3] | units.m, [3,2,1] | units.m)
        self.assertIsOfOrder([1,1,1] | units.m, 1 | units.m)
        self.assertIsOfOrder(1 | units.m, [0.4,1,3] | units.m)
        self.assertRaises(self.failureException, self.assertIsOfOrder, 
            1.0 | units.m, [0.3, 1.0, 4.0] | units.m, expected_message = 
            "@0, 1.0 m is not of order 0.3 m\n@2, 1.0 m is not of order 4.0 m")
    
    def test5(self):
        self.assertRaises(AmuseException, (1 | units.m).as_quantity_in, units.kg)
        self.assertRaises(AmuseException, (1 | units.m).as_quantity_in, 1 | units.cm)
        self.assertRaises(AmuseException, (1 | units.m).as_quantity_in, 1 | units.cm, 
            expected_message = "Cannot expres a unit in a quantity")
        # check that a failureException is raised, when we erroneously assert an exception:
        self.assertRaises(self.failureException, self.assertRaises, AmuseException, 
            (1 | units.m).as_quantity_in, units.cm, expected_message = 
            "AmuseException not raised")
        # the ultimate test... lol
        self.assertRaises(self.failureException, self.assertRaises, self.failureException, 
            self.assertRaises, AmuseException, (1 | units.m).as_quantity_in, 
            1 | units.cm, expected_message = "AssertionError not raised")
        self.assertRaises(AmuseException, lambda: 1 + (1|units.m))
        self.assertRaises(ZeroDivisionError, lambda: 1|units.m/0)
    
