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
       
