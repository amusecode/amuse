import unittest

from amuse.support.units.nbody import *

class TestNbodyUnits(unittest.TestCase):
   def test1(self):
       convert_nbody = nbody_to_si(units.parsec(1),units.MSun(20))
       y = mass(1)
       self.assertEqual(str(y), '1 nbody mass')
       y_in_si = convert_nbody.to_si(y)
       y_in_msun = y_in_si.in_(units.MSun)
       self.assertEqual(str(y_in_msun), '20.0 MSun')
       y_in_nbody = convert_nbody.to_nbody(y_in_msun) 
       self.assertEqual(str(y_in_nbody), '1.0 nbody mass')
       