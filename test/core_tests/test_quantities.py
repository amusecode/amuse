from amuse.test import amusetest

import numpy
import sys

from amuse.support.units import si, units, nbody_system
from amuse.support.data import core
from amuse.support.data.values import *

class TestQuantities(amusetest.TestCase):

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
        
    def test7(self):
        x = 2.0 | si.kg
        y = 1 / x
        self.assertEquals(y.value_in(1/si.kg), 0.5)
        
    
    def test8(self):
        x = (1.0, 2.0, 3.0) | si.kg
        self.assertTrue(isinstance(x, VectorQuantity))
        
    def test9(self):
        converter = nbody_system.nbody_to_si(1 | si.kg, 2 | si.s)
        self.assertEquals(0.0 | units.none, converter.to_nbody(zero))
        self.assertEquals(converter.to_nbody(zero), 0.0 | units.none)
    
    def test10(self):
        self.assertEquals(1 | units.m, 1 | units.m)
        self.assertTrue (1 | units.m == 1 | units.m)
        self.assertFalse(1 | units.m == 2 | units.m)
        self.assertTrue (1 | units.m != 2 | units.m)
        self.assertFalse(1 | units.m != 1 | units.m)
        self.assertTrue (1 | units.m >= 1 | units.m)
        self.assertFalse(1 | units.m >= 2 | units.m)
        self.assertTrue (1 | units.m <= 1 | units.m)
        self.assertFalse(1 | units.m <= 0 | units.m)
        self.assertTrue (1 | units.m >  0 | units.m)
        self.assertFalse(1 | units.m >  1 | units.m)
        self.assertTrue (1 | units.m <  3 | units.m)
        self.assertFalse(1 | units.m <  0 | units.m)
    
    def test11(self):
        self.assertEquals([1] | units.m, [1] | units.m)
        self.assertTrue ([1] | units.m == [1] | units.m)
        self.assertFalse([1] | units.m == [2] | units.m)
        self.assertTrue ([1] | units.m != [2] | units.m)
        self.assertFalse([1] | units.m != [1] | units.m)
        self.assertTrue ([1] | units.m >= [1] | units.m)
        self.assertFalse([1] | units.m >= [2] | units.m)
        self.assertTrue ([1] | units.m <= [1] | units.m)
        self.assertFalse([1] | units.m <= [0] | units.m)
        self.assertTrue ([1] | units.m >  [0] | units.m)
        self.assertFalse([1] | units.m >  [1] | units.m)
        self.assertTrue ([1] | units.m <  [3] | units.m)
        self.assertFalse([1] | units.m <  [0] | units.m)
    
    def test12(self):
        self.assertEquals(zero, zero)
        self.assertTrue (zero == zero)
        self.assertFalse(zero == zero + (1 | units.m))
        self.assertFalse(zero + (1 | units.m) == zero)
        self.assertTrue (zero != zero + (1 | units.m))
        self.assertFalse(zero != zero)
        self.assertTrue (zero >= zero)
        self.assertFalse(zero >= zero + (1 | units.m))
        self.assertTrue (zero <= zero + (1 | units.m))
        self.assertFalse(zero <= zero - (1 | units.m))
        self.assertTrue (zero >  zero - (1 | units.m))
        self.assertFalse(zero >  zero)
        self.assertTrue (zero <  zero + (1 | units.m))
        self.assertFalse(zero <  zero - (1 | units.m))
        self.assertTrue (zero == 0 | units.m)
    
    def test13(self):
        self.assertEquals('a' | units.string, 'a' | units.string)
        self.assertTrue ('a' | units.string == 'a' | units.string)
        self.assertFalse('a' | units.string == 'ab' | units.string)
        self.assertTrue ('a' | units.string != 'A' | units.string)
        self.assertFalse('a' | units.string != 'a' | units.string)
        self.assertTrue ('b' | units.string >= 'a' | units.string)
        self.assertFalse('B' | units.string >= 'a' | units.string)
        self.assertTrue ('a' | units.string <= 'ab' | units.string)
        self.assertFalse('a' | units.string <= 'A' | units.string)
        self.assertTrue ('a' | units.string >  'A' | units.string)
        self.assertFalse('a' | units.string >  'a' | units.string)
        self.assertTrue ('a' | units.string <  'b' | units.string)
        self.assertFalse('a' | units.string <  'B' | units.string)
       
        
        
