import unittest

from amuse.support.units.units import *

class TestUnitConversions(unittest.TestCase):
    def test1(self):
        km = 1000 * m
        self.assertEqual(1000, km.in_(m).number)
        self.assertEqual(0.001, m.in_(km).number)
    def test2(self):
        km = 1000 * m
        val = km(10)
        self.assertEqual(10000, val.in_(m).number)
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
        val = ms(10)
        self.assertEqual(36,val.in_(kmh).number)
    def test5(self):
        km = named('kilometer','km',1000 * m)
        h = named('hour','h',60 * 60 * s)
        kmh = km/h
        ms = m/s
        val = ms(10)
        self.assertEqual(36,val.in_(kmh).number)