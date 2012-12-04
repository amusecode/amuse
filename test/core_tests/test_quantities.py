from amuse.test import amusetest

import numpy
import sys


from amuse.support.exceptions import AmuseException
from amuse.units.quantities import *
from amuse.units import si
from amuse.units import units
from amuse.units import nbody_system
from amuse.units import quantities


from amuse import datamodel
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
        stars = datamodel.Particles(number_of_stars)
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
        self.assertEquals(0.0 | nbody_system.mass, converter.to_nbody(zero))
        self.assertEquals(converter.to_nbody(zero), 0.0 | nbody_system.mass)
    
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
        self.assertEquals('a', 'a')
        self.assertTrue ('a' == 'a')
        self.assertFalse('a' == 'ab')
        self.assertTrue ('a' != 'A')
        self.assertFalse('a' != 'a')
        self.assertTrue ('b' >= 'a')
        self.assertFalse('B' >= 'a')
        self.assertTrue ('a' <= 'ab')
        self.assertFalse('a' <= 'A')
        self.assertTrue ('a' >  'A')
        self.assertFalse('a' >  'a')
        self.assertTrue ('a' <  'b')
        self.assertFalse('a' <  'B')
    
    def test14(self):
        # Tests for 'is_quantity'
        self.assertTrue( is_quantity(0 | units.kg) )
        self.assertTrue( is_quantity(1 | units.none) )
        self.assertTrue( is_quantity([1.0, 2.0, 3.0] | units.m) )
        self.assertFalse( is_quantity(1) )
        self.assertFalse( is_quantity(1.0) )
        self.assertFalse( is_quantity("string") )
    
    def test15(self):
        # Tests for 'to_quantity'
        self.assertTrue( is_quantity(to_quantity(0 | units.kg)) )
        self.assertTrue( is_quantity(to_quantity(1 | units.none)) )
        self.assertTrue( is_quantity(to_quantity([1.0, 2.0, 3.0] | units.m)) )
        self.assertTrue( is_quantity(to_quantity(1)) )
        self.assertTrue( is_quantity(to_quantity(1.0)) )
        masses = [1, 2, 3] | units.kg
        self.assertTrue( to_quantity(masses) is masses )
        numbers = [1, 2, 3]
        self.assertFalse(          to_quantity(numbers) is numbers | units.none )
        self.assertTrue( numpy.all(to_quantity(numbers) == numbers | units.none) )
    
    def test16(self):
        # Tests for add/sub of quantity (with none unit) and number
        self.assertEquals( (2.0 | units.none) + 1.0,  3.0  )
        self.assertEquals( (2.0 | units.none) - 1.0,  1.0  )
        self.assertEquals( 1.0 + (2.0 | units.none),  3.0  )
        self.assertEquals( 1.0 - (2.0 | units.none), -1.0  )
    
    def test17(self):
        # Tests for add/sub of quantity (with other unit) and number
        number = 1.0
        quantity = 2.0 | units.m
        self.assertTrue( number.__add__(quantity) is NotImplemented)
        self.assertRaises(AmuseException, quantity.__radd__, number)
        self.assertTrue( number.__radd__(quantity) is NotImplemented)
        self.assertRaises(AmuseException, quantity.__add__, number)
        self.assertTrue( number.__sub__(quantity) is NotImplemented)
        self.assertRaises(AmuseException, quantity.__rsub__, number)
        self.assertTrue( number.__rsub__(quantity) is NotImplemented)
        self.assertRaises(AmuseException, quantity.__sub__, number)
        # in other words...
        self.assertRaises(AmuseException, lambda: number + quantity, 
            expected_message = "Cannot express none in m, the units do not have the same bases")
        self.assertRaises(AmuseException, lambda: quantity + number, 
            expected_message = "Cannot express none in m, the units do not have the same bases")
        self.assertRaises(AmuseException, lambda: number - quantity, 
            expected_message = "Cannot express none in m, the units do not have the same bases")
        self.assertRaises(AmuseException, lambda: quantity - number, 
            expected_message = "Cannot express none in m, the units do not have the same bases")
    
    def test18(self):
        quantity = 'string'
        self.assertEquals(quantity ,  'string')
        quantity = u'string'
        self.assertEquals(quantity , u'string')

    def test19(self):
        x = 1.0 | si.kg
        self.assertTrue(x==x.amin())
        self.assertTrue(x==x.prod())
        self.assertTrue(x==x.sorted())
        self.assertTrue(x==x.amax())
        self.assertTrue(x==x.sum())
    
    def test20(self):
        lengths = [] | units.m
        lengths.append(1 | units.m)
        self.assertEqual(lengths, [1] | units.m)
        lengths.append(2 | units.m)
        self.assertEqual(lengths, [1, 2] | units.m)
        
        positions = [] | units.m
        positions.append([1, 2, 3] | units.m)
        self.assertEqual(positions, [[1, 2, 3]] | units.m)
        positions.append([4, 5, 6] | units.m)
        self.assertEqual(positions, [[1, 2, 3], [4, 5, 6]] | units.m)
        
        two_positions = [] | units.m
        two_positions.append([[1, 2, 3], [-1, -2, -3]] | units.m)
        self.assertEqual(two_positions, [[[1, 2, 3], [-1, -2, -3]]] | units.m)
        two_positions.append([[4, 5, 6], [7, 8, 9]] | units.m)
        self.assertEqual(two_positions, [[[1, 2, 3], [-1, -2, -3]], [[4, 5, 6], [7, 8, 9]]] | units.m)
        # Appending quantities with incompatible shapes:
        two_positions.append(99 | units.m)
        self.assertEqual(two_positions, [1, 2, 3, -1, -2, -3, 4, 5, 6, 7, 8, 9, 99] | units.m)
    
    def test21(self):
        zero_vector = zero.as_vector_with_length(3)
        self.assertEqual(str(zero_vector), "[0.0, 0.0, 0.0] zero")
        
        self.assertEqual(zero_vector + (1 | units.m), [1, 1, 1] | units.m)
        self.assertEqual(zero_vector - (1 | units.m), [-1, -1, -1] | units.m)
        self.assertEqual((1 | units.m) + zero_vector, [1, 1, 1] | units.m)
        self.assertEqual((1 | units.m) - zero_vector, [1, 1, 1] | units.m)
        
        self.assertEqual(zero_vector + ([1, 1, 1] | units.m), [1, 1, 1] | units.m)
        self.assertEqual(zero_vector - ([1, 1, 1] | units.m), [-1, -1, -1] | units.m)
        self.assertEqual(([1, 1, 1] | units.m) + zero_vector, [1, 1, 1] | units.m)
        self.assertEqual(([1, 1, 1] | units.m) - zero_vector, [1, 1, 1] | units.m)
        
        for one_zero in zero_vector:
            self.assertEquals(one_zero, zero)
        
        self.assertEquals(zero_vector[0:2], zero.as_vector_with_length(2))
    

    def test22(self):
        x = numpy.asarray([1,2,3,4])
        y = 2 | units.m
        self.assertTrue(is_quantity(y * x))
        self.assertAlmostRelativeEquals(y*x, [2,4,6,8] | units.m)
        self.assertTrue(is_quantity(x * y))
        self.assertAlmostRelativeEquals(x*y, [2,4,6,8] | units.m)
    
    def test23(self):
        z = zero.as_vector_with_length(2)
        self.assertEquals(len(z), 2)
        z += 1 | units.kg
        self.assertEquals(z.unit, units.kg)
        self.assertEquals(z, [1,1] | units.kg)
    
    
    def xtest24(self):
        rhs = 2 | units.MSun / units.AU
        lhs = 3 | units.AU
        product = rhs * lhs
        product_unit = product.unit
        print product
        self.assertTrue(product_unit is units.MSun)
        
    
    def xtest25(self):
        rhs = 2.0 | (2 * units.MSun)**2 / units.AU
        lhs = 3.0 | units.MSun
        product = rhs / lhs
        product_unit = product.unit
        print product
        self.assertEquals(product_unit , units.MSun / units.AU)
        self.assertEquals(product_unit.local_factor , 4)
        
    
    def xtest26(self):
        rhs = 2.0 | units.AU / (2 * units.MSun)**2 
        lhs = 3.0 | units.MSun
        product = rhs * lhs
        product_unit = product.unit
        print product
        print type(product_unit)
        self.assertEquals(product_unit , units.AU / units.MSun)
        self.assertEquals(product_unit.local_factor , 1/4.0)
        
        
class TestAdaptingVectorQuantities(amusetest.TestCase):

    def test1(self):
        x = AdaptingVectorQuantity()
        self.assertEquals(x.append.__name__, "append_start")
        x.append(1 | units.kg)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 1)
        self.assertEquals(x.append.__name__, "append_normal")
        
        self.assertTrue(isinstance(x._number_list, list))
        self.assertFalse(isinstance(x._number_list, numpy.ndarray))
        self.assertTrue(isinstance(x.number, numpy.ndarray))
        self.assertFalse(isinstance(x.number, list))
        self.assertEquals(x._number_list, [1])
        self.assertEquals(x.number, numpy.array([1]))
    
    def test2(self):
        x = AdaptingVectorQuantity()
        self.assertEquals(len(x), 0)
        self.assertEquals(len(x.number), 0)
        self.assertEquals(str(x), '[]')
        x.append(1 | units.kg)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 1)
        self.assertEquals(str(x), '[1] kg')
        x.append(2 | units.kg)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 2)
        self.assertEquals(str(x), '[1, 2] kg')
    
    def test3(self):
        x = AdaptingVectorQuantity()
        x.extend([1,2,3] | units.kg)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 3)
        self.assertEquals(x.number, numpy.array([1,2,3]))
        
        x.extend([1,2,3] | units.g)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 6)
        self.assertAlmostRelativeEquals(x, [1000,2000,3000,1,2,3] | units.g)
    
    def test4(self):
        x = AdaptingVectorQuantity()
        x.prepend(1 | units.kg)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 1)
        self.assertEquals(str(x), '[1] kg')
        x.prepend(2 | units.kg)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 2)
        self.assertEquals(str(x), '[2, 1] kg')
    
    def test5(self):
        # Everything mixed...
        x = AdaptingVectorQuantity()
        x.extend([3,4,5] | units.kg)
        x.append(6 | units.kg)
        x.prepend(2 | units.kg)
        x.extend([7000,8000,9000] | units.g)
        x.prepend(1000 | units.g)
        x.append(10000 | units.g)
        self.assertEquals(x.unit, units.kg)
        self.assertEquals(len(x), 10)
        self.assertEquals(x.number, numpy.array([1,2,3,4,5,6,7,8,9,10]))
        self.assertEquals(x, [1,2,3,4,5,6,7,8,9,10]|units.kg)

class TestNumpyFunctionWithUnits(amusetest.TestCase):

    def test1(self):
        array = quantities.arange(0 | units.kg, 10 | units.kg, 1 | units.kg)
        self.assertEquals(len(array),  10)
        self.assertAlmostRelativeEquals(array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] | units.kg)
    
        
    def test2(self):
        array = quantities.linspace(0 | units.kg, 10 | units.kg, 11)
        self.assertEquals(len(array),  11)
        self.assertAlmostRelativeEquals(array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] | units.kg)
        
        
    def test3(self):
        array = quantities.linspace(0 , 10 , 11)
        self.assertEquals(len(array),  11)
        self.assertAlmostRelativeEquals(array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        
        
        
