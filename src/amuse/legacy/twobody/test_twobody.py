import amuse.legacy.twobody.twobody as twobody
from amuse.support.units import units
import unittest

class test_twobody(unittest.TestCase):

  def test1(self):
    nb=twobody.twobody()
    nb.new_particle(5.9742e24,6371.,0.,7.e6,-1.2124e7,0.,2.6679e3,4.6210e3)
    nb.evolve(3600.)
    state,err=nb.get_state(0)
    self.assertAlmostEqual(state['x'],0.,7)
    self.assertAlmostEqual(state['y']/(-3.30153815385e6),1.,7)
    self.assertAlmostEqual(state['z']/7.41119830598e6,1.,7)
    self.assertAlmostEqual(state['vx'],0.,7)
    self.assertAlmostEqual(state['vy']/(-8.29786882484e3),1.,7)
    self.assertAlmostEqual(state['vz']/(-0.967872571269e3),1.,7)
