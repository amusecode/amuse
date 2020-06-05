import numpy

from amuse.test import amusetest
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic import limepy

class TestLimepy(amusetest.TestCase):
    def test1(self):
          cluster=limepy.Limepy(5, 1, N=100).result

          self.assertAlmostEqual(cluster.total_mass().number, 1.0)
          self.assertEqual(len(cluster),100)

    def test2(self):
          cluster = limepy.new_limepy_model(7, 2, N=10000)

          self.assertAlmostEqual(cluster.total_mass().number, 1.0)
          self.assertEqual(len(cluster),10000)

    def test3(self):
          c = nbody_system.nbody_to_si(200 | units.MSun, 2 | units.parsec)
          cluster = limepy.new_limepy_model(7, 2, N=100, converter=c)

          self.assertAlmostEqual(cluster.total_mass(), 200. | units.MSun)
          self.assertEqual(len(cluster),100)

    def setUp(self):
          if not limepy.scipy_imported:
            self.skip("scipy not installed")
