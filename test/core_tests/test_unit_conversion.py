import unittest

from amuse.support.units.units import *

class TestUnitConversions(unittest.TestCase):
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
