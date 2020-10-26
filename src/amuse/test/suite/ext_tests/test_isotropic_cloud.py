from amuse.test import amusetest
from amuse.units import units, nbody_system, constants
from amuse.ic.isotropic_cloud import new_isotropic_cloud
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.datamodel import Particles

class TestSphericalIsotropicCloud(amusetest.TestCase):
  
  def test1(self):
    cloud = new_isotropic_cloud(2, m_cloud=0.|units.MSun)
    self.assertEqual(len(cloud),2)
    self.assertEqual(cloud.mass.sum(), 0.|units.MSun)
    
  def test2(self):
    m_star = 0.666|units.MSun
    a_min=666.|units.AU
    a_max=6666.|units.AU
    q_min=6.|units.AU
    cloud = new_isotropic_cloud(66,
                                m_star=m_star,
                                a_min=a_min,
                                a_max=a_max,
                                q_min=q_min)
    binary = Particles(1)
    binary[0].mass = m_star
    binary[0].position = (0.,0.,0.) | units.AU
    binary[0].velocity = (0.,0.,0.) | units.kms
    for comet in cloud:
      binary.add_particle(comet)
      mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per = \
        orbital_elements_from_binary(binary, G=constants.G)
      print(mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per)
      self.assertTrue( a_min < semimajor_axis < a_max )
      self.assertTrue( q_min < semimajor_axis*(1.-eccentricity) )
      binary.remove_particle(comet)
