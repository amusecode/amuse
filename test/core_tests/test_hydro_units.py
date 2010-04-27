from amuse.test import amusetest
from amuse.support.units.hydro_system import *
from amuse.support.units import constants

class TestHydroUnits(amusetest.TestCase):
    def test1(self):
        convert_hydro = hydro_to_si(1.0|units.g, 5.0|units.ms, 2.0|units.s)
        L = 1 | length
        L_in_si = convert_hydro.to_si(L)
        self.assertAlmostRelativeEqual(L_in_si, 10|units.m, 3)

    def test2(self):
        convert_hydro = hydro_to_si(10|units.N, 1|units.m, 1|units.s)
        F = 1 | force
        F_in_si = convert_hydro.to_si(F)
        self.assertAlmostRelativeEqual(F_in_si, 10|units.N, 3)

    def test4(self):
        convert_hydro = hydro_to_si(constants.G, 1|units.m, 1|units.s)
        M = 1 | mass
        M_in_si = convert_hydro.to_si(M)
        self.assertAlmostRelativeEqual(M_in_si, 1.49828895e+10|units.kg, 3)

    def test5(self):
        convert_hydro = hydro_to_si(constants.c, constants.hbar, constants.G)
        M = 1 | mass
        T = 1 | time
        L = 1 | length

        M_in_si = convert_hydro.to_si(M)
        T_in_si = convert_hydro.to_si(T)
        L_in_si = convert_hydro.to_si(L)

        self.assertAlmostRelativeEqual(M_in_si, 2.1764411e-8|units.kg, 3)
        self.assertAlmostRelativeEqual(T_in_si, 5.39124e-44|units.s, 3)
        self.assertAlmostRelativeEqual(L_in_si, 1.616252e-35|units.m, 3)


  
