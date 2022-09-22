import numpy

from amuse.test import amusetest
from amuse.units import units, nbody_system
from amuse.ic.brokenimf import *

# Instead of random, use evenly distributed numbers, just for testing
default_options = dict(random=False)

class TestMultiplePartIMF(amusetest.TestCase):
    
    def test1(self):
        print("Test MultiplePartIMF with default mass_boundaries and alphas, i.e. Salpeter")
        instance = MultiplePartIMF(mass_max=100.0 | units.MSun)
        self.assertEqual(instance.mass_boundaries, [0.1, 100.0] | units.MSun)
        self.assertEqual(instance.alphas, [-2.35])
        self.assertEqual(instance.number_of_bins, 1)
        self.assertEqual(instance.fraction_per_bin, [1.0])
        self.assertEqual(instance.cumulative_fractions, [0.0, 1.0])
        
        self.assertAlmostEqual(instance.mass([0.0]), 0.1 | units.MSun)
        self.assertAlmostEqual(instance.mass([1.0]), 100.0 | units.MSun)
        self.assertAlmostEqual(instance.mass_mean(), 0.351 | units.MSun, 3)
    
    def test2(self):
        print("Test MultiplePartIMF with mass_boundaries and alphas")
        instance = MultiplePartIMF(mass_boundaries = [1.0, 10.0, 100.0] | units.MSun, 
            alphas = [1.3, -3.3], **default_options)
        self.assertEqual(instance.mass_boundaries, [1.0, 10.0, 100.0] | units.MSun)
        self.assertEqual(instance.alphas, [1.3, -3.3])
        self.assertEqual(instance.number_of_bins, 2)
        self.assertAlmostEqual(instance.fraction_per_bin, numpy.array([0.5, 0.5]))
        self.assertEqual(instance.cumulative_fractions, [0.0, 0.5, 1.0])
        
        self.assertAlmostEqual(instance.mass([0.0]), 1.0 | units.MSun)
        self.assertAlmostEqual(instance.mass([0.5]), 10.0 | units.MSun)
        self.assertAlmostEqual(instance.mass([1.0]), 100.0 | units.MSun)
        
        self.assertAlmostEqual(instance.mass_mean(), 11.9457684987 | units.MSun)
        self.assertAlmostEqual(instance.mass_mean(), instance.next_mass(10000).mean(), 2)
    
    def test3(self):
        print("Test new_broken_power_law_mass_distribution with default mass_boundaries and alphas, i.e. Salpeter")
        masses = new_broken_power_law_mass_distribution(10000, mass_max=100.0 | units.MSun, **default_options)
        self.assertTrue((masses >= 0.1 | units.MSun).all())
        self.assertTrue((masses <= 100.0 | units.MSun).all())
        self.assertAlmostEqual(min(masses), 0.1 | units.MSun)
        self.assertAlmostEqual(max(masses), 100.0 | units.MSun)
        
        mass_mean = MultiplePartIMF(mass_boundaries=[0.1, 100.0]|units.MSun, 
            alphas=[-2.35]).mass_mean()
        self.assertAlmostEqual(mass_mean, 0.35136877959 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), 0.351 | units.MSun, 1)
    
    def test4(self):
        print("Test new_broken_power_law_mass_distribution with mass_boundaries and alphas")
        masses = new_broken_power_law_mass_distribution(10000, 
            mass_boundaries = [1.0, 10.0, 100.0] | units.MSun, 
            alphas = [1.3, -3.3], **default_options)
        self.assertTrue((masses >= 1.0 | units.MSun).all())
        self.assertTrue((masses <= 100.0 | units.MSun).all())
        self.assertAlmostEqual(min(masses), 1.0 | units.MSun)
        self.assertAlmostEqual(max(masses), 100.0 | units.MSun)
        
        mass_mean = MultiplePartIMF(mass_boundaries=[1.0, 10.0, 100.0]|units.MSun, 
            alphas=[1.3, -3.3]).mass_mean()
        self.assertAlmostEqual(mass_mean, 11.9457684987 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), 11.9457684987 | units.MSun, 1)
    
    def test5(self):
        print("Test new_scalo_mass_distribution")
        masses = new_scalo_mass_distribution(10000, **default_options)
        self.assertTrue((masses >= 0.1 | units.MSun).all())
        self.assertTrue((masses <= 125.0 | units.MSun).all())
        self.assertAlmostEqual(min(masses), 0.1 | units.MSun)
        self.assertAlmostEqual(max(masses), 125.0 | units.MSun)
        
        mass_mean = MultiplePartIMF(mass_boundaries=[0.10, 0.18, 0.42, 0.62, 1.18, 3.5, 125.0]|units.MSun, 
            alphas=[1.6, -1.01, -2.75, -2.08, -3.5, -2.63]).mass_mean()
        self.assertAlmostEqual(mass_mean, 0.487756751788 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), 0.487756751788 | units.MSun, 1)
    
    def test6(self):
        print("Test new_miller_scalo_mass_distribution")
        masses = new_miller_scalo_mass_distribution(10000, **default_options)
        self.assertTrue((masses >= 0.1 | units.MSun).all())
        self.assertTrue((masses <= 125.0 | units.MSun).all())
        self.assertAlmostEqual(min(masses), 0.1 | units.MSun)
        self.assertAlmostEqual(max(masses), 125.0 | units.MSun)
        
        mass_mean = MultiplePartIMF(mass_boundaries=[0.1, 1.0, 2.0, 10.0, 125.0]|units.MSun, 
            alphas=[-1.25, -2.0, -2.3, -3.3]).mass_mean()
        self.assertAlmostEqual(mass_mean, 0.885783055149 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), 0.885783055149 | units.MSun, 1)
    
    def test7(self):
        print("Test new_kroupa_mass_distribution")
        masses = new_kroupa_mass_distribution(10000, **default_options)
        self.assertTrue((masses >= 0.01 | units.MSun).all())
        roundoff = 1.0 + 1.0e-12
        self.assertTrue((masses <= (100.0 * roundoff) | units.MSun).all())
        self.assertAlmostEqual(min(masses), 0.01 | units.MSun)
        self.assertAlmostEqual(max(masses), 100.0 | units.MSun)
        
        mass_mean = MultiplePartIMF(mass_boundaries=[0.01, 0.08, 0.5, 100.0]|units.MSun, 
            alphas=[-0.3, -1.3, -2.3]).mass_mean()
        self.assertAlmostEqual(mass_mean, 0.376175542639 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), 0.376175542639 | units.MSun, 1)
    
    def test8(self):
        print("Test with problematic alphas (new_salpeter_mass_distribution would give zero division errors)")
        masses = new_broken_power_law_mass_distribution(10000, 
            mass_boundaries = [1.0, 10.0, 100.0] | units.MSun, 
            alphas = [-1, -2], **default_options)
        self.assertTrue((masses >= 1.0 | units.MSun).all())
        roundoff = 1.0 + 1.0e-12
        self.assertTrue((masses <= (100.0 * roundoff) | units.MSun).all())
        self.assertAlmostEqual(min(masses), 1.0 | units.MSun)
        self.assertAlmostEqual(max(masses), 100.0 | units.MSun)
        
        mass_mean = MultiplePartIMF(mass_boundaries=[1.0, 10.0, 100.0] | units.MSun, 
            alphas=[-1, -2]).mass_mean()
        self.assertAlmostEqual(mass_mean, 10.0 | units.MSun)
        self.assertAlmostRelativeEqual(masses.mean(), 10.0 | units.MSun, 1)
    
        masses = new_broken_power_law_mass_distribution(101, 
            mass_boundaries = [1.0, 100.0] | units.MSun, 
            alphas = [-1], **default_options)
        self.assertAlmostEqual(masses.median(), 10.0 | units.MSun)
    

