# simple cosmology calc.
# to be extended as needed

import numpy
from amuse.units import units, generic_unit_system
from amuse.units.quantities import to_quantity
from amuse.datamodel import Particles
from amuse.support.exceptions import AmuseException

def findbin(ylist,y):
  s=1
  if ylist[0]>=ylist[-1]:
    s=-1
  if s*y <= s*ylist[0]:
    return -1
  if s*y >= s*ylist[-1]:
    return len(ylist)    
  up=len(ylist)-1
  low=0
  while up-low>1:
    b=(low+up)/2
    if s*y < s*ylist[b]:
      up=b
    else:
      low=b
  return up       

class Hermitelookup(object):
  def __init__(self,xlist,ylist,yderiv):
    self.xlist=xlist
    self.ylist=ylist
    self.yderiv=yderiv

  def interpolatecubic(self,x, b):
    if b <= 0:
      return self.ylist[0]
    if b > len(self.ylist)-1:
      return self.ylist[-1]
    dx=self.xlist[b]-self.xlist[b-1]    
    dy=self.ylist[b]-self.ylist[b-1]    
    if dx==0.:
      return (self.ylist[b-1]+self.ylist[b])/2
    y1=self.ylist[b-1]
    yd2=self.yderiv[b]
    yd1=self.yderiv[b-1]
    u=(x-self.xlist[b-1])/(self.xlist[b]-self.xlist[b-1])
    return u**3*(-2*dy+dx*(yd1+yd2))+u**2*(3*dy-dx*(2*yd1+yd2))+dx*yd1*u+y1

  def evaluate(self,x):
    return self.interpolatecubic(x,findbin(self.xlist,x))

class Cosmology(object):
  def __init__(self,  # default=fifth year wmap+BAO+SN parameters, hinshaw 2008
                 omega=1.,
                 omegal = 0.726,
                 omegak = 0.,
                 omegar = 8.37e-5, # 4.165E-5/(h*h) includes 3 massless neutrino species, T0 = 2.72528
                 h = 0.705,
                 sigma8 = 0.812,
                 n=1000): 
    self.omega=omega
    self.omegal=omegal
    self.omegar=omegar
    self.omegak=omegak                   
    self.hubble0=h*(100 | units.kms/units.Mpc)
    self.omegam = omega - (omegak + omegar + omegal) 
    self.n=n
    
    a=(numpy.array(range(self.n+1))/float(self.n))**2
    t=[0.]
    dtda=[0.]
    dadt=[0.]
    for i in range(1,self.n+1):
      _t=t[-1]+1./6.*( self.invfriedmanint(a[i])+ 
                       self.invfriedmanint(a[i-1])+ 
                       4*self.invfriedmanint((a[i]+a[i-1])/2) )*(a[i]-a[i-1])
      t.append( _t )
      dtda.append(self.invfriedmanint(a[i]))
      dadt.append(1./dtda[-1])
    self.a=a
    self.t=numpy.array(t)
    self.dtda=numpy.array(dtda)
    self.dadt=numpy.array(dadt)
    self.age_lookup=Hermitelookup(self.a,self.t,self.dtda)
    self.a_lookup=Hermitelookup(self.t,self.a,self.dadt)
        
  def invfriedmanint(self,a):
    return a/(self.omegam*a+self.omegar+self.omegal*a**4+self.omegak*a**2)**0.5

  def hubble(self,a):
    return self.hubble0*self.dadtau(a)/a

  def dadtau(self,a):
    return (self.omegam/a+self.omegar/a**2+self.omegal*a**2+self.omegak)**0.5

  def d2dadtau2(self,a):
    return -1./2.*self.omegam/a**2-self.omegar/a**3+self.omegal*a
  
  def agefromz(self,z):
    return self.agefroma(1./(z+1.))

  def taufromz(self,z):
    return self.taufroma(1./(z+1.))
  
  def agefroma(self,a):
    return self.age_lookup.evaluate(a)/self.hubble0

  def taufroma(self,a):
    return self.age_lookup.evaluate(a)
    
  def afromage(self,age):
    return self.a_lookup.evaluate(age*self.hubble0)

  def afromtau(self,tau):
    return self.a_lookup.evaluate(tau)
    
def convert_comoving_to_physical(original, redshift=0.0, hubble_parameter=1.0, attribute_names=None):
    """
    Converts quantities or particle sets from comoving coordinates to physical 
    coordinates. In comoving coordinates, changes in positions due to the 
    expansion of the universe are corrected for, by dividing by the scale 
    factor 'a':
    a = 1 / (1 + z).
    This function will undo this correction for all quantities with units that 
    are (derived from) length units (e.g. position, velocity, energy, but not 
    time or mass).
    
    Optionally, a value for the Hubble parameter (value of Hubble constant in 
    units of 100 km/s/Mpc) can be supplied. If so, the units of 'original' are 
    assumed to be based on length/h, mass/h, and time/h, instead of length, 
    mass, and time, respectively. These factors will be divided out in the result
    """
    if isinstance(original, Particles):
        copy = original.copy_to_memory()
        if attribute_names is None:
            attribute_names = copy.get_attribute_names_defined_in_store()
        for attribute in attribute_names:
            setattr(copy, attribute, convert_quantity_from_comoving_to_physical(
                getattr(copy, attribute), redshift, hubble_parameter))
        return copy
    elif hasattr(original, "unit"):
        return convert_quantity_from_comoving_to_physical(original, redshift, hubble_parameter)
    else:
        raise AmuseException("Can't convert instance of {0} from comoving to physical "
            "coordinates (only Particles or Quantity supported)".format(original.__class__))

def convert_quantity_from_comoving_to_physical(original, redshift, hubble_parameter=1.0):
    for (exponent, unit) in to_quantity(original).unit.base:
        if unit is units.m or unit is generic_unit_system.length:
            return original * (1 / (1.0 + redshift))**exponent
    return original


if __name__=="__main__":  
  cosmo=Cosmology()
  print cosmo.agefroma(1.).in_(units.Myr)
  print cosmo.afromage(cosmo.agefroma(1.))


