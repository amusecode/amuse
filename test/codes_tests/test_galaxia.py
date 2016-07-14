from amuse.test import amusetest
from amuse.units import units, nbody_system, constants

from amuse.community.galaxia.interface import BarAndSpirals3D

class TestBarAndSpirals3D(amusetest.TestCase):
    def test1(self):
      galaxy= BarAndSpirals3D()
      self.assertEqual(galaxy.model_time, 0. | units.yr)
