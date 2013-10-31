import numpy

from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.units import nbody_system
from amuse.units import generic_unit_system
from amuse.units import units

from amuse import datamodel

class bb79_cloud(object):
  """
  spherical uniformly rotating cloud of particles with density perturbation (m=2)
  -- initial conditions for the 'standard isothermal test case' 
     Boss & Bodenheimer (1979, http://adsabs.harvard.edu/abs/1979ApJ...234..289B)
     -> binary fragmentation during isothermal collapse 

  arguments:
    targetN -- intended number of particles
    omega -- angular velocity (cloud rotates as a rigid body around the z-axis), 
      given in units of base rad/s if a converter is given, in 1./nbody_system.time if no converter is given
    rho_perturb -- amplitude of the density perturbation
    ethep_ratio -- ratio between total thermal and potential enegry
    convert_nbody -- to set the Nbody units
    base_grid -- base grid
  
  Default values set as in Boss & Bodenheimer (1979).
  
  In the first step, uniform sphere of particles is generated; then the azimuth of
  particles is changed (while radius is kept constant) to achieve the cosine density
  perturbation. See Kitsionas (2003, sec. 3.1, 
  http://adsabs.harvard.edu/abs/2003PhDT.......219K)
  """
  def __init__(self, targetN=10000, omega=0.775066020047 | nbody_system.time**-1, 
               rho_perturb=0.5, ethep_ratio=0.25, convert_nbody=None, base_grid=None):
    self.targetN=targetN
    if convert_nbody is not None:
      omega=convert_nbody.to_nbody(omega)
    self.omega=omega.value_in(1./nbody_system.time)
    self.rho_peturb=rho_perturb
    self.ethep_ratio=ethep_ratio
    self.convert_nbody=convert_nbody
    self.base_grid=base_grid
      
  def new_model(self):
        
    base_sphere=uniform_unit_sphere(self.targetN,base_grid=self.base_grid)
    x_uni,y_uni,z=base_sphere.make_xyz()
    self.actualN=len(x_uni)
    rad=numpy.sqrt(x_uni**2 + y_uni**2)
    phi=numpy.arctan2(y_uni,x_uni)
    n_vec=2000
    phi_new_vec=numpy.linspace(-numpy.pi, numpy.pi, n_vec)
    phi_old_vec=phi_new_vec + self.rho_peturb*(numpy.sin(2.*phi_new_vec)/2.)
    phi_new=numpy.interp(phi,phi_old_vec,phi_new_vec)
    x=rad*numpy.cos(phi_new)
    y=rad*numpy.sin(phi_new)
    
    rad=numpy.sqrt(x**2 + y**2)
    phi=numpy.arctan2(y,x)
    vel=self.omega*rad
    vx=-vel*numpy.sin(phi)
    vy= vel*numpy.cos(phi)
    vz=0.

    mass=numpy.ones_like(x)/self.actualN

    Ep=3./5
    self.internalE=Ep*self.ethep_ratio
    internal_energy=numpy.ones_like(x)*self.internalE
    
    return (mass,x,y,z,vx,vy,vz,internal_energy)

  @property
  def result(self):
    mass,x,y,z,vx,vy,vz,u = self.new_model()
    result = datamodel.Particles(self.actualN)
    result.mass = nbody_system.mass.new_quantity(mass)
    result.x = nbody_system.length.new_quantity(x)
    result.y = nbody_system.length.new_quantity(y)
    result.z = nbody_system.length.new_quantity(z)
    result.vx = nbody_system.speed.new_quantity(vx)
    result.vy = nbody_system.speed.new_quantity(vy)
    result.vz = nbody_system.speed.new_quantity(vz)
    result.u = (nbody_system.speed**2).new_quantity(u)

    if not self.convert_nbody is None:
        result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_generic())
        result = result.copy()

    return result
