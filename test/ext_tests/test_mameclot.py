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

