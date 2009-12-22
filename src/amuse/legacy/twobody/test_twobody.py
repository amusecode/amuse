import amuse.legacy.twobody.twobody as twobody
from amuse.support.units import units
import unittest
import numpy

class test_twobody(unittest.TestCase):

  def test_stumpff(self):
    self.assertAlmostEqual(twobody.stumpff_C(0),twobody.stumpff_C(0.0001),5)
    self.assertAlmostEqual(twobody.stumpff_C(0),twobody.stumpff_C(-0.0001),5)
    self.assertAlmostEqual(twobody.stumpff_S(0),twobody.stumpff_S(0.0001),5)
    self.assertAlmostEqual(twobody.stumpff_S(0),twobody.stumpff_S(-0.0001),5)

  def test1(self):
    nb=twobody.twobody()
    nb.new_particle(5.9742e24,6.371e6,0.,7.e6,-1.2124e7,0.,2.6679e3,4.6210e3)
    nb.evolve(3600.)
    state,err=nb.get_state(0)
    self.assertAlmostEqual(state['x'],0.,7)
    self.assertAlmostEqual(state['y']/(-3.30153815385e6),1.,7)
    self.assertAlmostEqual(state['z']/7.41119830598e6,1.,7)
    self.assertAlmostEqual(state['vx'],0.,7)
    self.assertAlmostEqual(state['vy']/(-8.29786882484e3),1.,7)
    self.assertAlmostEqual(state['vz']/(-0.967872571269e3),1.,7)

  def test2(self):
    nb=twobody.twobody()
    nb.new_particle(5.9742e24,7.1e6,0.,7.e6,-1.2124e7,0.,2.6679e3,4.6210e3)
    err=nb.evolve(3600.)
    self.assertEqual(err,1)
    dt,err=nb.get_time()
    self.assertAlmostEqual(dt/2584.9554627,1.,7)
    state,err=nb.get_state(0)
    self.assertAlmostEqual((state['x']**2+state['y']**2+state['z']**2)/(7.1e6)**2,1.,7)
