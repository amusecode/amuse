from amuse.test import amusetest

import numpy
import sys


from amuse.support.exceptions import AmuseException
from amuse.units.quantities import *
from amuse.units import si
from amuse.units import units
from amuse.units import trigo
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.units import core


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
        self.assertEqual(str(x[1]), "2.0 kg")


    def test4(self):
        g = si.kg / 1000
        x = [1.0, 2.0, 3.0] | si.kg
        self.assertEqual(str(x), "[1.0, 2.0, 3.0] kg")
        x[0] = 3000.0 | g
        self.assertEqual(str(x), "[3.0, 2.0, 3.0] kg")

    def test5(self):
        number_of_stars = 10
        stars = datamodel.Particles(number_of_stars)
        stars.position = [0,0,0] | units.km
        for i, star in enumerate(stars):
            star.position = units.km.new_quantity([float(i+1), float((i+1)*2), float(-1 * (i+1))])


        minpos = [float(sys.maxsize)] * 3 | units.m
        maxpos = [-float(sys.maxsize)] * 3 | units.m
        for star in stars:
            for i in range(3):
                if star.position[i] < minpos[i]:
                    minpos[i] = star.position[i]
                if star.position[i] > maxpos[i]:
                    maxpos[i] = star.position[i]

        self.assertEqual(str(minpos), "[1000.0, 2000.0, -10000.0] m")
        self.assertEqual(str(maxpos), "[10000.0, 20000.0, -1000.0] m")

    def test6(self):
        x = [1.0, 2.0, 3.0] | si.kg
        y = x.copy()
        y[0] = 3.0 | si.kg
        self.assertEqual(x[0].value_in(si.kg), 1.0)
        self.assertEqual(y[0].value_in(si.kg), 3.0)

    def test7(self):
        x = 2.0 | si.kg
        y = 1 / x
        self.assertEqual(y.value_in(1/si.kg), 0.5)


    def test8(self):
        x = (1.0, 2.0, 3.0) | si.kg
        self.assertTrue(isinstance(x, VectorQuantity))

    def test9(self):
        converter = nbody_system.nbody_to_si(1 | si.kg, 2 | si.s)
        self.assertEqual(0.0 | nbody_system.mass, converter.to_nbody(zero))
        self.assertEqual(converter.to_nbody(zero), 0.0 | nbody_system.mass)

    def test10(self):
        self.assertEqual(1 | units.m, 1 | units.m)
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
        self.assertEqual([1] | units.m, [1] | units.m)
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
        self.assertEqual(zero, zero)
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
        self.assertEqual('a', 'a')
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
        self.assertEqual( (2.0 | units.none) + 1.0,  3.0  )
        self.assertEqual( (2.0 | units.none) - 1.0,  1.0  )
        self.assertEqual( 1.0 + (2.0 | units.none),  3.0  )
        self.assertEqual( 1.0 - (2.0 | units.none), -1.0  )

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
        self.assertEqual(quantity ,  'string')
        quantity = 'string'
        self.assertEqual(quantity , 'string')

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
            self.assertEqual(one_zero, zero)

        self.assertEqual(zero_vector[0:2], zero.as_vector_with_length(2))


    def test22(self):
        x = numpy.asarray([1,2,3,4])
        y = 2 | units.m
        self.assertTrue(is_quantity(y * x))
        self.assertAlmostRelativeEquals(y*x, [2,4,6,8] | units.m)
        self.assertTrue(is_quantity(x * y))
        self.assertAlmostRelativeEquals(x*y, [2,4,6,8] | units.m)

    def test23(self):
        z = zero.as_vector_with_length(2)
        self.assertEqual(len(z), 2)
        z += 1 | units.kg
        self.assertEqual(z.unit, units.kg)
        self.assertEqual(z, [1,1] | units.kg)


    def xtest24(self):
        rhs = 2 | units.MSun / units.AU
        lhs = 3 | units.AU
        product = rhs * lhs
        product_unit = product.unit
        print(product)
        self.assertTrue(product_unit is units.MSun)


    def xtest25(self):
        rhs = 2.0 | (2 * units.MSun)**2 / units.AU
        lhs = 3.0 | units.MSun
        product = rhs / lhs
        product_unit = product.unit
        print(product)
        self.assertEqual(product_unit , units.MSun / units.AU)
        self.assertEqual(product_unit.local_factor , 4)


    def xtest26(self):
        rhs = 2.0 | units.AU / (2 * units.MSun)**2
        lhs = 3.0 | units.MSun
        product = rhs * lhs
        product_unit = product.unit
        print(product)
        print(type(product_unit))
        self.assertEqual(product_unit , units.AU / units.MSun)
        self.assertEqual(product_unit.local_factor , 1/4.0)

    def test27(self):
        a=[1.|units.kg,2.|units.kg,3000.| units.g, 4.| (1000*units.g)]
        b=VectorQuantity.new_from_scalar_quantities(*a)
        c=[1.,2.,3.,4.] | units.kg
        print(a[0].unit==a[2].unit)
        self.assertEqual(b,c)

    def test28(self):
        a=[1.|units.kg,2.|units.kg,3000.| units.m, 4.| (1000*units.g)]
        try:
          b=VectorQuantity.new_from_scalar_quantities(*a)
          raise Exception("expect error")
        except:
          pass
    
    def test29(self):
        one_inch = 2.54 | units.cm
        self.assertFalse(isinstance(one_inch, core.unit))
        self.assertTrue(isinstance(one_inch.as_unit(), core.unit))
        self.assertEqual(one_inch.as_unit(), 2.54 * units.cm)
        self.assertEqual(1 | one_inch.as_unit(), 2.54 | units.cm)
    
    def test30(self):
        a=1.5| units.km
        b=1000. | units.m
        self.assertEqual(a%b,500. | units.m)
        a=[1.5,1.75]| units.km
        b=1000. | units.m
        self.assertEqual(a%b,[500.,750] | units.m)
        a=[1.5,1.75]| units.km
        b=[1000.,500.] | units.m
        self.assertEqual(a%b,[500.,250.] | units.m)
        a=[1.5]| units.km
        b=[1000.,500.] | units.m
        self.assertEqual(a%b,[500.,0.] | units.m)

    def test31(self):
        """ 
        test trigonometric unit stuff
        """
        self.assertEqual(units.pi,numpy.pi)
        a=units.pi
        self.assertEqual(trigo.to_rad(a), numpy.pi | units.rad)
        self.assertEqual(trigo.to_deg(a), 180. | units.deg)
        self.assertEqual(trigo.to_rev(a), 0.5 | units.rev)
        a=90 | units.deg
        self.assertEqual(trigo.to_rad(a), numpy.pi/2 | units.rad)
        self.assertEqual(trigo.to_deg(a), 90. | units.deg)
        self.assertEqual(trigo.to_rev(a), 0.25 | units.rev)
        a=0.75 | units.rev
        self.assertEqual(trigo.to_rad(a), 3/2.*numpy.pi | units.rad)
        self.assertEqual(trigo.to_deg(a), 270. | units.deg)
        self.assertEqual(trigo.to_rev(a), 0.75 | units.rev)
        a=2*numpy.pi
        self.assertEqual(trigo.to_rad(a), 2*numpy.pi | units.rad)
        self.assertEqual(trigo.to_deg(a), 360. | units.deg)
        self.assertEqual(trigo.to_rev(a), 1. | units.rev)

        a=45. | units.deg
        self.assertEqual(trigo.sin(a),numpy.sin(45./180*numpy.pi))
        self.assertEqual(trigo.cos(a),numpy.cos(45./180*numpy.pi))
        self.assertEqual(trigo.tan(a),numpy.tan(45./180*numpy.pi))

        a=1. | units.rad
        self.assertEqual(trigo.sin(a),numpy.sin(1.))
        self.assertEqual(trigo.cos(a),numpy.cos(1.))
        self.assertEqual(trigo.tan(a),numpy.tan(1.))

        a=0.125 | units.rev
        self.assertEqual(trigo.sin(a),numpy.sin(45./180*numpy.pi))
        self.assertEqual(trigo.cos(a),numpy.cos(45./180*numpy.pi))
        self.assertEqual(trigo.tan(a),numpy.tan(45./180*numpy.pi))

        a=45. | units.deg
        self.assertAlmostEqual(trigo.arcsin(trigo.sin(a)),45. | units.deg,13)
        self.assertAlmostEqual(trigo.arccos(trigo.cos(a)),45. | units.deg,13)
        self.assertAlmostEqual(trigo.arctan(trigo.tan(a)),45. | units.deg,13)

    def test32(self):
        a=numpy.array([[1.,2.,3.],[4.,5.,6.]]) |  units.m
        b=numpy.array([[1.,2.,3.],[4.,5.,6.]])
        self.assertEqual(list(a.flatten()), list(a.flat))
        flat1=b.flat
        flat2=a.flat
        self.assertEqual(flat1[2:5],flat2[2:5].number)
        next(flat1)
        next(flat2)
        self.assertEqual(flat1.index,flat2.index)
        self.assertEqual(flat1.base,flat2.base.number)
        self.assertEqual(flat1.copy(),flat2.copy().number)        

    def test32b(self):
        a=numpy.array([[1.,2.,3.],[4.,5.,6.]]) |  units.m
        b=numpy.array([[1.,2.,3.],[4.,5.,6.]])
        flat1=b.flat
        flat2=a.flat
        self.assertEqual(flat1[2:5],flat2[2:5].number)
        self.assertEqual(flat1,flat2.number)
        
        flat2[:]=numpy.arange(6) | units.cm
        a_=numpy.array([[0.,1.,2.],[3.,4.,5.]]) |  units.cm
        self.assertEqual(a, a_)

    def test33(self):
        a=[1,2,3,4]
        b=new_quantity_nonone(a,units.none)
        self.assertEqual(len(a),len(b))
        b=new_quantity_nonone(a,2*units.none)
        self.assertEqual(len(a),len(b))

    def test34(self):
        a= [1,2,3,4,5] | units.m
        x= a.value_in(units.cm)
        x[0]=-1
        self.assertEqual(a, [1,2,3,4,5] | units.m)

        a= [1,2,3,4,5] | units.m
        x= a.value_in(units.m)
        x[0]=-1
        self.assertEqual(a, [1,2,3,4,5] | units.m)


class TestAdaptingVectorQuantities(amusetest.TestCase):

    def test1(self):
        x = AdaptingVectorQuantity()
        self.assertEqual(x.append.__name__, "append_start")
        x.append(1 | units.kg)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 1)
        self.assertEqual(x.append.__name__, "append_normal")

        self.assertTrue(isinstance(x._number_list, list))
        self.assertFalse(isinstance(x._number_list, numpy.ndarray))
        self.assertTrue(isinstance(x.number, numpy.ndarray))
        self.assertFalse(isinstance(x.number, list))
        self.assertEqual(x._number_list, [1])
        self.assertEqual(x.number, numpy.array([1]))

    def test2(self):
        x = AdaptingVectorQuantity()
        self.assertEqual(len(x), 0)
        self.assertEqual(len(x.number), 0)
        self.assertEqual(str(x), '[]')
        x.append(1 | units.kg)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 1)
        self.assertEqual(str(x), '[1] kg')
        x.append(2 | units.kg)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 2)
        self.assertEqual(str(x), '[1, 2] kg')

    def test3(self):
        x = AdaptingVectorQuantity()
        x.extend([1,2,3] | units.kg)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 3)
        self.assertEqual(x.number, numpy.array([1,2,3]))

        x.extend([1,2,3] | units.g)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 6)
        self.assertAlmostRelativeEquals(x, [1000,2000,3000,1,2,3] | units.g)

    def test4(self):
        x = AdaptingVectorQuantity()
        x.prepend(1 | units.kg)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 1)
        self.assertEqual(str(x), '[1] kg')
        x.prepend(2 | units.kg)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 2)
        self.assertEqual(str(x), '[2, 1] kg')

    def test5(self):
        # Everything mixed...
        x = AdaptingVectorQuantity()
        x.extend([3,4,5] | units.kg)
        x.append(6 | units.kg)
        x.prepend(2 | units.kg)
        x.extend([7000,8000,9000] | units.g)
        x.prepend(1000 | units.g)
        x.append(10000 | units.g)
        self.assertEqual(x.unit, units.kg)
        self.assertEqual(len(x), 10)
        self.assertEqual(x.number, numpy.array([1,2,3,4,5,6,7,8,9,10]))
        self.assertEqual(x, [1,2,3,4,5,6,7,8,9,10]|units.kg)

    def test6(self):
        x =  6 | units.kg
        y =  5 | units.kg
        self.assertEqual(x/y, 6/5) 
        self.assertEqual(x//y, 6//5) 
        self.assertEqual(operator.__truediv__(x,y), 1.2)

class TestNumpyFunctionWithUnits(amusetest.TestCase):

    def test1(self):
        array = quantities.arange(0 | units.kg, 10 | units.kg, 1 | units.kg)
        self.assertEqual(len(array),  10)
        self.assertAlmostRelativeEquals(array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] | units.kg)


    def test2(self):
        array = quantities.linspace(0 | units.kg, 10 | units.kg, 11)
        self.assertEqual(len(array),  11)
        self.assertAlmostRelativeEquals(array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] | units.kg)


    def test3(self):
        array = quantities.linspace(0 , 10 , 11)
        self.assertEqual(len(array),  11)
        self.assertAlmostRelativeEquals(array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    def test4(self):
        x = quantities.arange(0 | units.yr, 10 | units.yr, 1 | units.yr)
        y = (2.0|units.km) * (x/ (2.0|units.yr))**2 + (20.0|units.km)

        fit = quantities.polyfit(x, y, 2)

        self.assertEqual(len(fit), 3)
        self.assertEqual(fit[0].unit, units.km/units.yr**2)

        fit_values = quantities.polyval(fit, x)

        self.assertEqual(fit_values.shape, x.shape)
        self.assertEqual(y.unit, fit_values.unit)

        self.assertAlmostRelativeEquals(y, fit_values, 1)

    def test5(self):
        a=[1,2,3] | units.m
        b=[4,5,6] | units.m
        
        ab1=quantities.column_stack((a,b))
        ab2=quantities.column_stack((a.number,b.number)) | units.m
        
        self.assertEqual(ab1,ab2)
        

