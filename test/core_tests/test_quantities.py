import unittest

from amuse.support.units import si
from amuse.support.data.values import *

class TestQuantities(unittest.TestCase):

    def test1(self):
        x = 1.0 | si.kg
        self.assertTrue(isinstance(x, ScalarQuantity))
        x = [1.0, 2.0, 3.0] | si.kg
        self.assertTrue(isinstance(x, VectorQuantity))
        
    def test2(self):
        x = [1.0, 2.0, 3.0] | si.kg
        y = [2.0, 3.0, 4.0] | si.kg
        xy = x * y
        self.assertTrue(isinstance(xy, VectorQuantity))