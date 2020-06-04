import numpy

from amuse.test import amusetest
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.mameclot import mameclot

class TestPlummer(amusetest.TestCase):
    def test1(self):
          cluster=mameclot().result
    
          self.assertAlmostEqual(cluster.total_mass().number, 1.0)
          self.assertEqual(len(cluster),10000)
          
    def test2(self):
          clusters,c1,c2=mameclot().result_split
    
          self.assertAlmostEqual((c1.total_mass()+c2.total_mass()).number, 1.0)
          self.assertEqual(len(c1),10000)          
          self.assertEqual(len(c2),0)

    def test3(self):
          clusters,c1,c2=mameclot(mass_ratio=0.25).result_split
    
          self.assertAlmostEqual((c1.total_mass()+c2.total_mass()).number, 1.0)
          self.assertEqual(len(c1),8000)          
          self.assertEqual(len(c2),2000)                    

