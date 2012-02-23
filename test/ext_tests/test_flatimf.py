import numpy

from amuse.test import amusetest
from amuse.units import units, nbody_system

from amuse.ic.flatimf import FlatIMF
from amuse.ic.flatimf import new_flat_mass_distribution, new_flat_mass_distribution_nbody

class FlatIMFTests(amusetest.TestCase):
    
    def test1(self):
        instance = FlatIMF(0.1 | units.MSun, 100.0 | units.MSun)
        self.assertAlmostEqual(instance.mass_mean(), 99.9/numpy.log(1e3) | units.MSun)
        instance = FlatIMF(42.0 | units.MSun, 42.00001 | units.MSun)
        self.assertAlmostEqual(instance.mass_mean(), 42.000 | units.MSun, 3)
    
    def test2(self):
        instance = FlatIMF(0.1 | units.MSun, 1000.0 | units.MSun)
        self.assertAlmostEqual(instance.mass(1.0), 1000.0 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.75), 100.0 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.5),   10.0 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.25),   1.0 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.0),    0.1 | units.MSun)
    
    def test3(self):
        numpy.random.seed(345672)
        instance = FlatIMF(0.1 | units.MSun, 100 | units.MSun)
        n = 10000
        total_mass, set_of_masses = instance.next_set(n)
        
        mean_mass = 99.9/numpy.log(1e3) | units.MSun
        self.assertAlmostEqual(instance.mass_mean(), mean_mass)
        self.assertAlmostRelativeEqual(total_mass / n, mean_mass, 3)
        self.assertAlmostEqual(total_mass / n, 14.4615334306 | units.MSun)
        self.assertAlmostEqual(instance.mass(0.5),   set_of_masses.median(), 1)
    
    def test4(self):
        print "Test 4: testing user interface (SI units)"
        numpy.random.seed(345672)
        masses = new_flat_mass_distribution(10000)
        
        self.assertEqual(len(masses), 10000)
        self.assertAlmostRelativeEqual(masses.mean(), FlatIMF().mass_mean(), 3)
        self.assertAlmostEqual(masses.mean(), 17.5145247111 | units.MSun)
        self.assertAlmostEqual(masses.amin(), 0.100145673289 | units.MSun)
        self.assertAlmostEqual(masses.amax(), 124.94980234 | units.MSun)
    
    def test5(self):
        print "Test 5: testing user interface (SI units), optional args"
        numpy.random.seed(345672)
        masses = new_flat_mass_distribution(10000, 
            mass_min=10.0|units.MSun, mass_max=100.0|units.MSun)
        
        self.assertEqual(len(masses), 10000)
        self.assertAlmostRelativeEqual(masses.mean(), 
            FlatIMF(mass_min=10.0|units.MSun, mass_max=100.0|units.MSun).mass_mean(), 2)
        self.assertAlmostEqual(masses.mean(), 39.1111546565 | units.MSun)
        self.assertAlmostEqual(masses.amin(), 10.0047015091 | units.MSun)
        self.assertAlmostEqual(masses.amax(), 99.9870310764 | units.MSun)
    
    def test6(self):
        print "Test 6: testing user interface (nbody units)"
        numpy.random.seed(345672)
        n = 10000
        masses = new_flat_mass_distribution_nbody(n)
        
        self.assertEqual(len(masses), n)
        self.assertAlmostEqual(masses.sum(),      1.0 | nbody_system.mass)
        self.assertAlmostEqual(masses.mean(), 1.0 / n | nbody_system.mass)
        self.assertAlmostEqual(masses.amin(), 0.100145673289 / 17.5145247111 / n | nbody_system.mass)
        self.assertAlmostEqual(masses.amax(), 124.94980234 / 17.5145247111 / n | nbody_system.mass)
    
    def test7(self):
        print "Test 7: testing user interface (nbody units), optional args"
        numpy.random.seed(345672)
        n = 10000
        masses = new_flat_mass_distribution_nbody(n, 
            mass_min=10.0|nbody_system.mass, mass_max=100.0|nbody_system.mass)
        
        self.assertEqual(len(masses), n)
        self.assertAlmostEqual(masses.sum(),      1.0 | nbody_system.mass)
        self.assertAlmostEqual(masses.mean(), 1.0 / n | nbody_system.mass)
        self.assertAlmostEqual(masses.amin(), 10.0047015091 / 39.1111546565 / n | nbody_system.mass)
        self.assertAlmostEqual(masses.amax(), 99.9870310764 / 39.1111546565 / n | nbody_system.mass)
    

