import unittest

from amuse.support.units import units
from amuse.support.data import core

class TestBase(unittest.TestCase):
    pass
   
class TestTemporalAttribute(TestBase):
    def test1(self):
        x = core.TemporalAttribute('mass')
        x.set_value_at_time(1 | units.s, 1 | units.m)
        x.set_value_at_time(2 | units.s, 2 | units.m)
        x.set_value_at_time(2.5 | units.s, 3 | units.m)
        x.set_value_at_time(3 | units.s, 4 | units.m)
        x.set_value_at_time(4 | units.s, 5 | units.m)
        x.set_value_at_time(5 | units.s, 6 | units.m)
        x.set_value_at_time(7 | units.s, 7 | units.m)
        self.assertEquals(x.get_value_at_time(4 | units.s)[1], 5 | units.m)
        self.assertEquals(x.get_value_at_time(3 | units.s)[1], 4 | units.m)
        self.assertEquals(x.get_value_at_time(2.5 | units.s)[1], 3 | units.m)
        self.assertEquals(x.get_value_at_time(7 | units.s)[1], 7 | units.m)
        self.assertEquals(x.get_value_at_time(1 | units.s)[1], 1 | units.m)
        
    def test1(self):
        x = core.TemporalAttribute('mass')
        for i in range(1, 2001):
            x.set_value_at_time(i | units.s, i | units.m)
        self.assertEquals(x.get_value_at_time(325 | units.s)[1], 325 | units.m)
        self.assertEquals(x.get_value_at_time(2000 | units.s)[1], 2000 | units.m)
        self.assertEquals(x.get_value_at_time(2010 | units.s)[1], 2000 | units.m)
