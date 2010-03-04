import sys
import os

from amuse.test import amusetest
from amuse.support.units import units
from amuse.ext.salpeter import SalpeterIMF
        
class SalpeterIMFTests(amusetest.TestCase):
    def test1(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass_mean().value_in(units.MSun), 0.351, 3)

    def test2(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass(1.0).value_in(units.MSun), 100, 3)
        self.assertAlmostEqual(instance.mass(0).value_in(units.MSun), 0.1, 3)
       
    def test3(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        n = 10000
        total_mass, set_of_masses = instance.next_set(10000)
        mean = total_mass.value_in(units.MSun) / float(n)
        exact_mean = instance.mass_mean().value_in(units.MSun)
        self.assertTrue(abs(mean - exact_mean) < 0.1)
        
    def test4(self):
        instance = SalpeterIMF(0.1 | units.MSun, 125 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual( 1.0 / instance.mass_mean().value_in(units.MSun), 2.8253, 4)
       

