from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import constants

from amuse.test import amusetest
import os.path

class Test(amusetest.TestCase):
    
    def Ntest1(self):
        I = constants.GetConstantsFromFiles()
        #I.get_table_from_url()
        I.get_table_from_file()
        I.check_current_file_with_table()
        
    def Ntest2(self):
        I = constants.Constants()
        I.generate_nist()
        
    def test3(self):
        from amuse.support.units import nist
        self.assertAlmostEqual(nist.h.value_in(nist.J*nist.s), 6.6e-34, 35)
        self.assertAlmostEqual(nist.c.value_in(nist.m/nist.s), 299792458.0, 7)
        
