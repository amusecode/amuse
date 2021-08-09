from amuse.test import amusetest

from amuse.units import units
from amuse.ic.brokenimf import new_kroupa_mass_distribution


class TestsForIssue777(amusetest.TestCase):
    def test_upper_segment(self):
        "Test if a star in the upper mass segment will get the right mass"
        lower_limit = 1.0 | units.MSun
        upper_limit = 1.0 | units.MSun
        mass = new_kroupa_mass_distribution(
            1,
            mass_min=lower_limit,
            mass_max=upper_limit,
        )
        self.assertEqual(mass[0], 1.0 | units.MSun)

    def test_middle_segment(self):
        "Test if a star in the middle mass segment will get the right mass"
        lower_limit = 0.2 | units.MSun
        upper_limit = 0.2 | units.MSun
        mass = new_kroupa_mass_distribution(
            1,
            mass_min=lower_limit,
            mass_max=upper_limit,
        )
        self.assertEqual(mass[0], 0.2 | units.MSun)

    def test_lower_segment(self):
        "Test if a star in the lower mass segment will get the right mass"
        lower_limit = 0.02 | units.MSun
        upper_limit = 0.02 | units.MSun
        mass = new_kroupa_mass_distribution(
            1,
            mass_min=lower_limit,
            mass_max=upper_limit,
        )
        self.assertEqual(mass[0], 0.02 | units.MSun)
