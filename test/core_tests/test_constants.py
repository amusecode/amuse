from amuse.units import units
from amuse.units import nist
from amuse.test import amusetest

import os.path

from amuse.support import data
class Test(amusetest.TestCase):
    
    def generate_constants_if_necessary(self):
        try:
            from amuse.support.units import constants
        except ImportError:
            I = nist.Constants()
            I.generate_constants()
            print "constants.py has been generated."
    
    def Ntest1(self):
        I = nist.GetConstantsFromFiles()
        #I.get_table_from_url()
        I.get_table_from_file()
        I.check_current_file_with_table()
        
    def Ntest2(self):
        I = nist.Constants()
        I.generate_constants()
        
    def test3(self):
        self.generate_constants_if_necessary()
        from amuse.support.units import constants
        self.assertAlmostEqual(constants.h.value_in(units.J*units.s), 6.6e-34, 35)
        self.assertAlmostEqual(constants.c.value_in(units.m/units.s), 299792458.0, 7)

    def test4(self):
        self.generate_constants_if_necessary()
        from amuse.support.units import constants
        self.assertAlmostEqual(constants.h, 6.6e-34 | units.J*units.s, 35)
        self.assertAlmostEqual(constants.c, 299792458.0 | units.m/units.s, 7)

    def test5(self):
        self.generate_constants_if_necessary()
        from amuse.support.units import constants
        self.assertAlmostEqual(constants.Planck_length**2 / 
            (constants.hbar*constants.G/constants.c**3), 1.0 | units.none, 5)
        self.assertAlmostEqual(constants.Planck_mass**2 / 
            (constants.hbar*constants.c/constants.G), 1.0 | units.none, 5)
        self.assertAlmostEqual(constants.Planck_time**2 / 
            (constants.hbar*constants.G/constants.c**5), 1.0 | units.none, 5)
        

