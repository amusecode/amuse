from amuse.test import amusetest
from amuse.units import units
from amuse.units.constants import *


class TestConstants(amusetest.TestCase):
    
    def test1(self):
        self.assertAlmostEqual(h.value_in(units.J*units.s), 6.6e-34, 35)
        self.assertAlmostEqual(c.value_in(units.m/units.s), 299792458.0, 7)

    def test2(self):
        self.assertAlmostEqual(h, 6.6e-34 | units.J*units.s, 35)
        self.assertAlmostEqual(c, 299792458.0 | units.m/units.s, 7)

    def test3(self):
        self.assertAlmostRelativeEquals(
            Planck_length**2,
            hbar * G / c**3, 5)
        self.assertAlmostRelativeEquals(
            Planck_mass**2,
            hbar * c / G, 5)
        self.assertAlmostRelativeEquals(
            Planck_time**2,
            hbar * G / c**5, 5)
        
    def test4(self):
        self.assertAlmostRelativeEquals(Rydberg_constant_times_hc_in_J, 1 | units.Ry, 7)
        self.assertAlmostRelativeEquals(2 * h * Rydberg_constant, 
            fine_hyphen_structure_constant**2 * electron_mass * c, 7)
        self.assertAlmostRelativeEquals(fine_hyphen_structure_constant, 
            elementary_charge**2 / (2 * h * c * electric_constant), 7)
    

