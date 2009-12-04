import unittest

import numpy
import sys

from amuse.support.units import si, units
from amuse.support.data import core
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
        
        
    def test3(self):
        x = [1.0, 2.0, 3.0] | si.kg
        y = [2.0, 3.0, 4.0] | si.kg
        self.assertTrue(isinstance(x[0], ScalarQuantity))
        self.assertEquals(str(x[1]), "2.0 kg")
        
        
    def test4(self):
        g = si.kg / 1000
        x = [1.0, 2.0, 3.0] | si.kg
        self.assertEquals(str(x), "[1.0, 2.0, 3.0] kg")
        x[0] = 3000.0 | g
        self.assertEquals(str(x), "[3.0, 2.0, 3.0] kg")
        
    def test5(self):
        number_of_stars = 10
        stars = core.Stars(number_of_stars)
        stars.position = [0,0,0] | units.km
        for i, star in enumerate(stars):
            star.position = units.km.new_quantity([float(i+1), float((i+1)*2), float(-1 * (i+1))])
        
        
        minpos = [float(sys.maxint)] * 3 | units.m
        maxpos = [-float(sys.maxint)] * 3 | units.m
        for star in stars:
            for i in range(3):
                print star.position[i], minpos[i]
                if star.position[i] < minpos[i]:
                    minpos[i] = star.position[i]
                if star.position[i] > maxpos[i]:
                    maxpos[i] = star.position[i]
    
        self.assertEquals(str(minpos), "[1000.0, 2000.0, -10000.0] m")
        self.assertEquals(str(maxpos), "[10000.0, 20000.0, -1000.0] m")
        
    def test6(self):
        x = [1.0, 2.0, 3.0] | si.kg
        y = x.copy()
        y[0] = 3.0 | si.kg
        self.assertEquals(x[0].value_in(si.kg), 1.0)
        self.assertEquals(y[0].value_in(si.kg), 3.0)
        
        
