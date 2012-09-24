from amuse.test import amusetest
import numpy

from amuse.support.exceptions import AmuseException



from amuse.units import core
from amuse.units.units import *
from amuse.units.constants import *
class TestUnitConversions(amusetest.TestCase):
    def test1(self):
        km = 1000 * m
        self.assertEqual(1000, km.value_in(m))
        self.assertEqual(0.001, m.value_in(km))
        
    def test2(self):
        km = 1000 * m
        val = km(10)
        self.assertEqual(10000, val.value_in(m))
        
    def test3(self):
        km = 1000 * m
        val = km(10)
        self.assertEqual("10 1000 * m", str(val))
        
        km = named('kilometer','km',1000 * m)
        val = km(10)
        self.assertEqual("10 km", str(val))
        
    def test4(self):
        km = named('kilometer','km',1000 * m)
        h = named('hour','h',60 * 60 * s)
        kmh = km/h
        ms = m/s
        val = 10 | m/s
        self.assertEqual(36,val.value_in(kmh))
        
    def test5(self):
        km = named('kilometer','km',1000 * m)
        h = named('hour','h',60 * 60 * s)
        kmh = km/h
        ms = m/s
        val = 10 | ms
        self.assertEqual(36,val.value_in(kmh))
        
    def test6(self):
        no1 = m / m
        no2 = no_unit
        self.assertEqual("m / m", str(no1))
        self.assertTrue(no1.has_same_base_as(no2))

        
    def test7(self):
        x = (100 * kg ** 2)  / kg 
        self.assertEqual("100 * kg**2 / kg", str(x))
        self.assertEqual("100.0 * kg", str(x.to_simple_form()))
        
        
    def test8(self):
        x = (10 | g) * (2 | m) 
        self.assertEqual("20 0.001 * m * kg", str(x))
        self.assertEqual("0.02 kg * m", str(x.as_quantity_in(kg * m)))
        x = (10 | kg) * (2000 | g) 
        self.assertEqual("20000 0.001 * kg**2", str(x))
        
    def test9(self):
        speed_of_light = 1 | (lightyear * yr**-1)
        time = 1e-9 | s
        length = speed_of_light * time
        length_in_m = length.value_in(m)
        self.assertAlmostEqual(0.2997988, length_in_m, 6)
        
    def test10(self):
        eps0_1 = mu0**-1*c**-2
        eps0_2 = (1 | none)/(mu0*(c**2))
        self.assertTrue((eps0_1 - eps0_2) / eps0_1 < (1e-10 | none))
        b =((1.|e)**2)
        f = (hbar*c*4.*numpy.pi* eps0)**-1
        fine_structure_constant_calculated = (b * f)
        fine_structure_constant = 7.297352537650e-3
        self.assertAlmostEquals(fine_structure_constant_calculated, fine_structure_constant, 5)
        
    

    def test11(self):
        vel1 = 1 | m / s
        self.assertRaises(core.IncompatibleUnitsException, vel1.as_quantity_in, s / m, 
            expected_message = "Cannot express m / s in s / m, the units do not have the same bases")
    
    def test12(self):
        self.assertEquals((1234 | g).as_string_in(g), '1234.0 g')
        self.assertEquals((1234 | g).as_string_in(kg), '1.234 kg')
        self.assertEquals((1.234 | kg).as_string_in(g), '1234.0 g')
        self.assertEquals((1.0 | km * s**-1).as_string_in(m / s), '1000.0 m / s')
        self.assertEquals((1.0 | km * s**-1).as_string_in(s**-1 * m), '1000.0 s**-1 * m')
        self.assertEquals((1.0 | km / s).as_string_in((10*J/g)**0.5), '10.0 (10 * J / g)**0.5')

class TestNonNumericUnits(amusetest.TestCase):
    def test1(self):
        string1 = "string" | string
        self.assertRaises(AmuseException, string1.as_quantity_in, m, 
            expected_message = "Cannot convert non-numeric quantities in to another unit")

    def test2(self):
        x = "test" | string
        self.assertEquals("test", x.value_in(string))  
            
    def test3(self):
        test_unit = core.enumeration_unit(
            "test", 
            "test",
            [1,2,3],
            ["one", "two", "three"]
        )
        x = 1 | test_unit
        self.assertEquals(1, x.value_in(test_unit))    
        self.assertEquals("one", str(x))  
        self.assertRaises(Exception, lambda: 4 | test_unit, 
            expected_message = "<4> is not a valid value for unit<test>")
    
    def test4(self):
        self.assertRaises(Exception, lambda: 1 | string, 
            expected_message = "<1> is not a valid value for unit<string>")
        
    def test5(self):
        test_unit = core.enumeration_unit(
            "test", 
            "test",
            [1,2,3],
            ["one", "two", "three"]
        )
        self.assertEquals(3, len(list(test_unit.quantities())))  
        for x, y in  zip(test_unit.quantities(), ["one", "two", "three"]):
            self.assertEquals(str(x), y)
    
       
    def test6(self):
        test_unit = core.enumeration_unit(
            "test", 
            "test",
            [1,4,7]
        )
        self.assertEquals(3, len(list(test_unit.quantities())))  
        for x, y in  zip(test_unit.quantities(), ["1", "4", "7"]):
            self.assertEquals(str(x), y)
    
    def test7(self):
        test_unit = core.enumeration_unit(
            "test", 
            "test",
            range(5)
        )
        self.assertEquals(5, len(list(test_unit.quantities())))
        self.assertEquals(1 | test_unit, 1 | test_unit)
        self.assertTrue (1 | test_unit == 1 | test_unit)
        self.assertFalse(1 | test_unit == 2 | test_unit)
        self.assertTrue (1 | test_unit != 2 | test_unit)
        self.assertFalse(1 | test_unit != 1 | test_unit)
        self.assertTrue (1 | test_unit >= 1 | test_unit)
        self.assertFalse(1 | test_unit >= 2 | test_unit)
        self.assertTrue (1 | test_unit <= 1 | test_unit)
        self.assertFalse(1 | test_unit <= 0 | test_unit)
        self.assertTrue (1 | test_unit >  0 | test_unit)
        self.assertFalse(1 | test_unit >  1 | test_unit)
        self.assertTrue (1 | test_unit <  3 | test_unit)
        self.assertFalse(1 | test_unit <  0 | test_unit)
        
    
