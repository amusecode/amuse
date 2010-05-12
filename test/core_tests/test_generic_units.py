from amuse.test import amusetest
from amuse.support.units.generic_unit_system import *
from amuse.support.units.generic_unit_converter import *
from amuse.support.units import constants

class TestGenericUnits(amusetest.TestCase):
    def test1(self):
        L = 1 | length
        T = 1 | time
        V = 1 | speed
        self.assertTrue(L/T == V)

    def test2(self):
        #natural units
        convert_generic = generic_to_si(constants.c, constants.hbar, constants.G, constants.kB)
        M = 1 | mass
        T = 1 | time
        L = 1 | length
        H = 1 | temperature

        M_in_si = convert_generic.to_si(M)
        T_in_si = convert_generic.to_si(T)
        L_in_si = convert_generic.to_si(L)
        H_in_si = convert_generic.to_si(H)

        self.assertAlmostRelativeEqual(M_in_si, 2.1764411e-8|units.kg, 3)
        self.assertAlmostRelativeEqual(T_in_si, 5.39124e-44|units.s, 3)
        self.assertAlmostRelativeEqual(L_in_si, 1.616252e-35|units.m, 3)
        self.assertAlmostRelativeEqual(H_in_si, 1.416785e32|units.K, 3)
