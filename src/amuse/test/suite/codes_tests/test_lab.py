import os
import sys
import numpy

from numpy import random
from amuse.test import amusetest
from amuse.test.amusetest import get_path_to_results


from amuse import lab
class TestCodesInLab(amusetest.TestCase):

    
    def test1(self):
        codes = (
            'SeBa',
            'BHTree',
            'Mocassin',
            'SPHRay',
            'Fi',
            'Gadget2',
            'Bonsai',
            'Octgrav',
        )
        
        for code in codes:
            self.assertTrue(hasattr(lab, code), msg = 'excpected class {0} in lab module, but the class was not found'.format(code))

