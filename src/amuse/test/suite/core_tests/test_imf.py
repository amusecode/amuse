import numpy.random

from amuse.test import amusetest
from amuse.units import units

from amuse.ic.brokenimf import new_kroupa_mass_distribution


class TestMassFunctions(amusetest.TestCase):

    def test_kroupa_default_masslimits(self):
        """
        Test the default mass limits of a Kroupa mass distribution
        """
        numpy.random.seed(124)
        mass_min = 0.01 | units.MSun
        mass_max = 125 | units.MSun
        mass = new_kroupa_mass_distribution(10000)
        self.assertTrue(mass.max() <= mass_max)
        self.assertTrue(mass.min() >= mass_min)

    def test_kroupa_lower_masslimit(self):
        """
        Test the lower mass limits of a Kroupa mass distribution
        """
        numpy.random.seed(124)
        mass_min = 1.0 | units.MSun
        mass = new_kroupa_mass_distribution(
            10000, mass_min=mass_min)
        self.assertTrue(mass.min() >= mass_min)

    def test_kroupa_upper_masslimit(self):
        """
        Test the upper mass limits of a Kroupa mass distribution
        """
        numpy.random.seed(124)

        mass_max = 1.0 | units.MSun
        mass = new_kroupa_mass_distribution(
            10000, mass_max=mass_max)
        self.assertTrue(mass.max() <= mass_max)

    def test_kroupa_masslimits(self):
        """
        Test both mass limits of a Kroupa mass distribution
        """
        numpy.random.seed(124)

        mass_max = 1.0 | units.MSun
        mass_min = 1.0 | units.MSun
        mass = new_kroupa_mass_distribution(
            10000, mass_min=mass_min, mass_max=mass_max)
        self.assertEqual(mass.max(), mass.min())
