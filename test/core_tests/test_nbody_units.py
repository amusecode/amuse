import unittest

from amuse.support.units.nbody_system import *

class TestNbodyUnits(unittest.TestCase):
    def test1(self):
       convert_nbody = nbody_to_si(units.parsec(1),units.MSun(20))
       y = 1 | mass
       self.assertEqual(str(y), '1 nbody mass')
       y_in_si = convert_nbody.to_si(y)
       y_in_msun = y_in_si.in_(units.MSun)
       self.assertEqual(str(y_in_msun), '20.0 MSun')
       y_in_nbody = convert_nbody.to_nbody(y_in_msun) 
       self.assertEqual(str(y_in_nbody), '1.0 nbody mass')
       
    def test2(self):
       convert_nbody = nbody_to_si(units.MSun(1.0), units.km(149.5e6))
       y = 29800 | units.m / units.s
       y_in_nbody = convert_nbody.to_nbody(y) 
       self.assertEqual(str(y_in_nbody.unit), 'nbody length * nbody time**-1')
       self.assertAlmostEqual(y_in_nbody.number, 1.0, 3)
       
    def test3(self):
       convert_nbody = nbody_to_si(units.MSun(1.0), units.km(149.5e6))
       y = 1 | length * (time**-1)
       y_in_si = convert_nbody.to_si(y) 
       #self.assertEqual(str(y_in_nbody.unit), 'nbody length * nbody time**-1')
       self.assertAlmostEqual(y_in_si.number, 29795.4, 1)
       
    def test4(self):
       convert_nbody = nbody_to_si(units.MSun(1.0), units.km(149.5e6))
       y_in_nbody = convert_nbody.to_nbody(units.G) 
       self.assertEqual(str(y_in_nbody.unit), 'nbody length**3 * nbody mass**-1 * nbody time**-2' )
       self.assertAlmostEqual(y_in_nbody.number, 1.0, 9)
       
    def test5(self):
       convert_nbody = nbody_to_si(units.km(149.5e6), units.MSun(1.0))
       y = 29800 | units.m / units.s
       y_in_nbody = convert_nbody.to_nbody(y) 
       self.assertEqual(str(y_in_nbody.unit), 'nbody length * nbody time**-1')
       self.assertAlmostEqual(y_in_nbody.number, 1.0, 3)
       

    def test6(self):
       convert_nbody = nbody_to_si(10 | units.kg , 5 | units.m / units. s )
       y = 5 | units.m / units.s
       y_in_nbody = convert_nbody.to_nbody(y) 
       self.assertEqual(str(y_in_nbody.unit), 'nbody length * nbody time**-1')
       self.assertAlmostEqual(y_in_nbody.number, 1.0, 3)
       y_in_si = convert_nbody.to_si(y_in_nbody) 
       self.assertAlmostEqual(y_in_si.number, 5.0, 3)
       
    def test7(self):
       convert_nbody = nbody_to_si(1 | units.kg , 1 | units.m / units. s )
       y = 1 | time
       y_in_si = convert_nbody.to_si(y) 
       self.assertEqual(str(y_in_si.unit), 's')
       self.assertAlmostEqual(y_in_si.number, 6.6730000000000296e-11, 3)
