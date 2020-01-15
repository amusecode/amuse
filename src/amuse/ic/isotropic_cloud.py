"""
Spherical isotropic (Oort) cloud -- initial conditions for a particle set

See for example:
  * Duncan, M.; Quinn, T.; Tremaine, S. -- http://adsabs.harvard.edu/abs/1987AJ.....94.1330D
  * Feng, F.; Bailer-Jones, C. A. L. -- http://adsabs.harvard.edu/abs/2014MNRAS.442.3653F
  * Dybczynski, P. A. -- http://adsabs.harvard.edu/abs/2002A%26A...396..283D

Initial conditions -- given by distributions of orbital elements:
  semi-major axes -- power law, default: dn/da ~ a^(-1.5),
  eccentricities -- dn/de ~ e,
  constrain on the minimum pericenter,
  isotropic orbits -- distribution of orbital inclinations: cos(i) = -1--1,
      longitude of ascending node: 0--2pi,
      argument of periastron: 0--2pi,
      mean anomaly: 0--2pi,
  equal mass particles.
"""

import numpy

from amuse.units import units, nbody_system
from amuse.datamodel import Particles
from amuse.community.kepler.interface import Kepler

def random_power_min_max(size, x_min, x_max, exp_plus_one):
  """
  returns random floats in the interval [x_min,x_max] drawn from distribution
  pdf(x) = const * x**(exp_plus_one-1), x_min <= x <= x_max; 
  assuming: x_min < x_max, exp_plus_one != 0
  """
  r = numpy.random.random(size=size)
  x_min_gamma = x_min**exp_plus_one
  x_max_gamma = x_max**exp_plus_one
  return (x_min_gamma + (x_max_gamma - x_min_gamma)*r)**(1./exp_plus_one)

def relative_position_and_velocity_from_orbital_elements(mass1,
                                                         mass2,
                                                         semimajor_axis,
                                                         eccentricity,
                                                         mean_anomaly,
                                                         seed=None):
  """
  Function that returns relative positions and velocity vectors or orbiters with masses 
  mass2 of the central body with mass mass1 in Cartesian coordinates;
  for vectors of orbital elements -- semi-major axes, eccentricities, mean anomalies.
  3D orientation of orbits (inclination, longitude of ascending node and argument of periapsis) are random.
  (cos(incl) is uniform -1--1, longitude of ascending node and argument of periapsis are uniform 0--2pi)
  Assuming mass1 is static in the center [0,0,0] m, [0,0,0] km/s (that is mass2<<mass1)
  """
  position_vectors = []
  velocity_vectors = []
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  kepler = Kepler(converter)
  kepler.initialize_code()
  r_vec = (0.,0.,0.) | units.AU
  v_vec = (0.,0.,0.) | units.kms
  # to change seed for each particle
  if seed is not None:
    i=0
  for m2_i, a_i, ecc_i, ma_i in zip(mass2, semimajor_axis, eccentricity, mean_anomaly):
    #print m2_i, a_i, ecc_i, ma_i
    if seed is not None:
      kepler.set_random(seed+i)
      i=i+1
    kepler.initialize_from_elements(mass=(mass1+m2_i),semi=a_i,ecc=ecc_i,mean_anomaly=ma_i,random_orientation=-1)
    ri = kepler.get_separation_vector()
    vi = kepler.get_velocity_vector()
    # this is to get ~half of the orbits retrograde (that is with inclination
    # of 90--180 degrees) --> velocity = -velocity
    vel_vec_dir = numpy.random.random()
    if (vel_vec_dir<=0.5):
      vel_orientation = 1.
    else:
      vel_orientation = -1.
    position_vectors.append([ri[0], ri[1], ri[2]])
    velocity_vectors.append([vel_orientation*vi[0], vel_orientation*vi[1], vel_orientation*vi[2]])
  kepler.stop()
  return position_vectors, velocity_vectors

def ecc_random_power_with_min_peri(n, semi, min_peri, power=2.):
  """
  random distribution in eccentricity P(e)~e
  (power = actual_exponent + 1)
  with minimum pericenter of min_peri
  for given semi-major axes semi
  """
  x = numpy.random.power(power,size=n)
  peri = semi*(1.-x)
  while numpy.any(peri<min_peri):
    filter_small_peri = (peri<min_peri)
    n_new = sum(filter_small_peri)
    #print "\t updating q", peri.min().in_(units.AU), n_new
    x_random_new = numpy.random.power(power,size=n_new)
    x_new = 1.*x
    x_new[filter_small_peri] = x_random_new
    x = 1.*x_new
    peri = semi*(1.-x)
  return x


class SphericalIsotropicCloud(object):
  def __init__(self,
               targetN,
               m_star=1.|units.MSun,
               m_cloud=0.|units.MSun,
               a_min=3000.|units.AU,
               a_max=10000.|units.AU,
               q_min=32.|units.AU,
               gamma=-1.5,
               seed=None):
    self.targetN = targetN
    self.m_star = m_star
    self.m_cloud = m_cloud
    self.a_min = a_min
    self.a_max = a_max
    self.q_min = q_min
    self.gamma = gamma
    self.seed = seed
    
    if (self.q_min is not None):
      if (self.q_min > self.a_min):
        self.a_min_q_corr = self.q_min
      else:
        self.a_min_q_corr = self.a_min
    else:
      self.a_min_q_corr = self.a_min
  
  def new_model(self):
    if self.seed is not None:
      numpy.random.seed(self.seed)
    
    a_in_au = random_power_min_max(self.targetN, 
                                   self.a_min_q_corr.value_in(units.AU),
                                   self.a_max.value_in(units.AU),
                                   self.gamma+1.)
    a = a_in_au * 1.|units.AU
    ecc = ecc_random_power_with_min_peri(self.targetN, a, self.q_min, power=2.)
    mean_anomaly = 2.*numpy.pi * numpy.random.random(size=self.targetN)
    m_comets = (self.m_cloud / self.targetN) * numpy.ones_like(ecc)
    position_vectors, velocity_vectors = \
      relative_position_and_velocity_from_orbital_elements(self.m_star,
                                                           m_comets,
                                                           a,
                                                           ecc,
                                                           mean_anomaly,
                                                           seed=self.seed)
    return (m_comets, position_vectors, velocity_vectors)
  
  @property
  def result(self):
    masses, position_vectors, velocity_vectors = self.new_model()
    result = Particles(self.targetN)
    result.mass = masses
    result.position = position_vectors
    result.velocity = velocity_vectors
    return result
  
def new_isotropic_cloud(number_of_particles, *list_arguments, **keyword_arguments):
  """
  Spherical isotropic cloud ~ Oort cloud given by distributions of orbital elements:
    semi-major axes -- power law, default: dn/da ~ a^(-1.5),
    eccentricities -- dn/de ~ e,
    constrain on the minimum pericenter,
    isotropic orbits -- distribution of orbital inclinations: cos(i) = -1--1,
        longitude of ascending node: 0--2pi,
        argument of periastron: 0--2pi,
        mean anomaly: 0--2pi,
    equal mass particles
  
  The default values correspond to papers:
  * Duncan, M.; Quinn, T.; Tremaine, S. -- http://adsabs.harvard.edu/abs/1987AJ.....94.1330D
  * Feng, F.; Bailer-Jones, C. A. L. -- http://adsabs.harvard.edu/abs/2014MNRAS.442.3653F (their DQT model)
  
  :argument number_of_particles: number of particles to include in the cloud
  :argument m_star: mass of the central star
  :argument m_cloud: total mass of the cloud (particles are equal mass)
  :argument a_min: minimal semimajor axis
  :argument a_max: maximal semimajor axis, a_min < a_max
  :argument q_min: minimal pericenter
  :argument gamma: exponent of the semimajor axis distribution, f(a) ~ a^(gamma)
  :argument seed: random seed -- set to reproduce exactly the same IC
  """
  uc = SphericalIsotropicCloud(number_of_particles, *list_arguments, **keyword_arguments)
  return uc.result

if __name__ in ('__main__'):
  cloud = new_isotropic_cloud(10)
  print(cloud)
