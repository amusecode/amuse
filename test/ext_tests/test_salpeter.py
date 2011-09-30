import sys
import os
import numpy.random

from amuse.test import amusetest
from amuse.units import units, nbody_system
from amuse.ic.salpeter import _SalpeterIMF
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

class SalpeterIMFTests(amusetest.TestCase):
    
    def test1(self):
        instance = _SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass_mean().value_in(units.MSun), 0.351, 3)
    
    def test2(self):
        instance = _SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass(1.0), 100 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.0), 0.1 | units.MSun)
    
    def test3(self):
        numpy.random.seed(345672)
        instance = _SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        n = 10000
        total_mass, set_of_masses = instance.next_set(n)
        
        self.assertAlmostEqual(instance.mass_mean(), 0.35136877959 | units.MSun)
        self.assertAlmostEqual(total_mass / n,       0.35136877959 | units.MSun, 1)
        self.assertAlmostEqual(total_mass / n,       0.33999456911 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.5),   set_of_masses.median(), 2)
    
    def test4(self):
        instance = _SalpeterIMF(0.1 | units.MSun, 125 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass_mean(), 0.353943475903 | units.MSun)
    
    def test5(self):
        print "Test 5: testing user interface (SI units)"
        numpy.random.seed(345672)
        masses = new_salpeter_mass_distribution(1000)
        
        self.assertEqual(len(masses), 1000)
        self.assertAlmostEqual(masses.mean(), 0.334475937397 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), _SalpeterIMF().mass_mean(), 1)
        self.assertAlmostEqual(masses.amin(), 0.10017909529 | units.MSun)
        self.assertAlmostEqual(masses.amax(), 19.7132849297 | units.MSun)
    
    def test6(self):
        print "Test 6: testing user interface (nbody units)"
        numpy.random.seed(345672)
        masses = new_salpeter_mass_distribution_nbody(1000)
        
        self.assertEqual(len(masses), 1000)
        self.assertAlmostEqual(masses.sum(),      1.0 | nbody_system.mass)
        self.assertAlmostEqual(masses.mean(), 1.0 / 1000 | nbody_system.mass)
        self.assertAlmostEqual(masses.amin(), 0.10017909529 / 0.334475937397 / 1000 | nbody_system.mass)
        self.assertAlmostEqual(masses.amax(), 19.7132849297 / 0.334475937397 / 1000 | nbody_system.mass)
    

