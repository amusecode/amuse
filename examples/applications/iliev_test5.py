"""
Radiative transfer comparison test 5 (Iliev et al 2009) 

This example demonstrates the coupling of a radiative transfer code to a hydrodynamics code to do
radiative hydrodynamic problems. Specifically it calculates the Iliev radiative trasnfer 
comparison project test 5 (expanding HII region in homogeneous medium)

It is simplified somewhat in the sense that it uses the SimpleX without any multifrequency 
radiative transfer or thermal evolution. This is not a fundamental limitation (the method would 
work with full thermal evolution available).

The coupling interface is with simple functions.

currently fi and simplex are hardcoded as base classes for SPH and radiative transfer

"""

import numpy
from amuse.support.units import units
from amuse.support.units import nbody_system

from amuse.support.units import constants
from amuse.community.simplex.interface import SimpleX
from amuse.community.fi.interface import Fi

from amuse.support.data import core

from amuse.ext.evrard_test import uniform_random_unit_cube,uniform_unit_sphere

try:
    from amuse import plot
    from matplotlib import pyplot
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False

from amuse.support.io import write_set_to_file


"""
 these define the properties of the medium in front and behind the ionisation front
 (necessary because no thermal evolution)
"""
mu=1.| units.amu
muion =0.5 | units.amu
xtrans=0.06
mutrans=mu/(1+xtrans)
Tinit=100. | units.K
Ttrans=13500. | units.K 
Tion=13500. | units.K
rhoinit=0.001 | (units.amu / units.cm**3)
uinit=constants.kB * Tinit/mu
utrans=constants.kB * Ttrans/mutrans
uion=constants.kB * Tion/muion
mp=None

def fake_internal_energy_from_xion(xion):
  """
  this function fakes an internal energy from the ionisation  
  (this is the main shortcut here)
  """
#  return uinit+(uion-uinit)*xion 
  u=uinit+(utrans-uinit)*xion/xtrans
  a=numpy.where(xion > xtrans )[0]
  u[a]=utrans+(uion-utrans)*(xion[a]-xtrans)/(1.-xtrans)
  return u

  
def glass(N, target_rms=0.05):
  """
   make glass for initial condition generation
  """  
  if target_rms < 0.001:
    print "warning: target_rms highly unlikely to succeed"
    
  L=1| nbody_system.length
  dt=0.01 | nbody_system.time
  x,y,z=uniform_random_unit_cube(N).make_xyz()
  vx,vy,vz=uniform_unit_sphere(N).make_xyz()
   

  p=core.Particles(N)
  p.x=L*x
  p.y=L*y
  p.z=L*z
  p.h_smooth=0. * L
  p.vx= 0.1*vx | (nbody_system.speed)
  p.vy= 0.1*vy | (nbody_system.speed)
  p.vz= 0.1*vz | (nbody_system.speed)
  p.u= (0.1*0.1) | nbody_system.speed**2 
  p.mass=(8./N) | nbody_system.mass

  sph=Fi(use_gl=False,mode='periodic',redirection='none')   
  sph.initialize_code()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.gamma=1
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=dt  
  sph.parameters.verbosity=0
  sph.parameters.pboxsize=2*L
  sph.parameters.artificial_viscosity_alpha=1.| units.none
  sph.parameters.beta=2. | units.none
  sph.commit_parameters()
  sph.gas_particles.add_particles(p)
  sph.commit_particles()

#  sph.start_viewer()

  t=0. | nbody_system.time
  rms=1.
  i=0
  while rms > target_rms:
    t=t+(0.25 | nbody_system.time)
    sph.evolve_model(t)
    h=sph.particles.h_smooth.value_in(nbody_system.length)
    rho=h**(-3.)
    rms=rho.std()/rho.mean()
    print rms

  x=sph.particles.x.value_in(nbody_system.length)
  y=sph.particles.y.value_in(nbody_system.length)
  z=sph.particles.z.value_in(nbody_system.length)

  del sph  
  return x,y,z


def iliev_test_5_ic( N=10000,
                  Ns=10,
                  L=15. | units.kpc ):
  """
  iliev test 5 particle distributions
   N= number of gas part
   Ns= number of sources (recommended to be order 10 for smoothrad distrbiution
   L=half boxsize
  """

  mp=rhoinit*(2*L)**3/N
 
#  x,y,z=uniform_random_unit_cube(N).make_xyz()
  x,y,z=glass(N,target_rms=0.05)
  
  p=core.Particles(N)
  p.x=L*x
  p.y=L*y
  p.z=L*z
  p.h_smooth=0. | units.parsec
  p.vx= 0. | (units.km/units.s)
  p.vy= 0. | (units.km/units.s)
  p.vz= 0. | (units.km/units.s)
  p.u= uinit 
  p.rho = rhoinit
  p.mass=mp
  p.flux=0. | (units.s**-1)
  p.xion=0. | units.none

  sources=core.Particles(Ns)
  x,y,z=uniform_unit_sphere(Ns).make_xyz()

  sources.x=L*x*(1./N)**(1./3)/10
  sources.y=L*y*(1./N)**(1./3)/10
  sources.z=L*z*(1./N)**(1./3)/10
  sources.rho=rhoinit/100.
  sources.flux=(5.e48/Ns) | (units.s**-1)
  sources.xion=1. | units.none
  
  return p,sources

iliev_test_5_ic.__test__ = False

def iliev_test_5( N=10000,
                  Ns=10,
                  L=15. | units.kpc,
                  dt=None):
  """
  prepare iliev test and return SPH and simplex interfaces
  """  
  gas,sources=iliev_test_5_ic(N,Ns,L)                


  conv=nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
     
  sph=Fi(conv,use_gl=False,mode='periodic',redirection='none')   
  sph.initialize_code()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.gamma=1
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=dt  
  sph.parameters.verbosity=0
  sph.parameters.pboxsize=2*L
  sph.commit_parameters()
  sph.gas_particles.add_particles(gas)
  sph.commit_particles()

#  sph.start_viewer()
         
  rad=SimpleX(number_of_workers = 1,redirection='none')
  rad.initialize_code()

  rad.parameters.box_size=2*L
  rad.parameters.hilbert_order=0

  rad.commit_parameters()

  gas.add_particles(sources)
  rad.particles.add_particles(gas)
  rad.commit_particles()
                  
  return sph,rad                  

iliev_test_5.__test__ = False

def update_pos_rho(sys,pa):
  p=pa.copy()
  p.rho=3./4/numpy.pi*8*p.mass/p.radius**3
  channel = p.new_channel_to(sys.particles)
  channel.copy_attributes(["x","y","z","rho"])

def update_u(sys,pa):
  p=pa.copy()
  p.u=fake_internal_energy_from_xion(p.xion)
  channel = p.new_channel_to(sys.particles)
  channel.copy_attributes(["u"])
  
def radhydro_evolve(sph,rad,tend,dt):
  """
   evolve function to co-evolve sph under radiative feedback from rad

   the evolve proceeds as follows (a form of leapfrog integrator):

   - 1/2 step sph
   - update positions and densities in rad
   - full step rad transfer
   - update internal energies in sph
   - 1/2 step sph 

   tend=end time
   dt=time step 

   this function dump snapshots in file dump-..

  """
  i=0
  write_set_to_file(sph.gas_particles,"dump-%6.6i"%i,"amuse",
                        append_to_file=False)    

  t=0. | units.Myr
  while t<tend-dt/2:
    print t        
    print "sph1"
    sph.evolve_model(t+dt/2)    
    print "rad"
    update_pos_rho(rad,sph.gas_particles)
    rad.evolve_model(t+dt)
    update_u(sph,rad.particles)
    print "sph2"
    sph.evolve_model(t+dt)    
    t+=dt
    i+=1
    write_set_to_file(sph.gas_particles,"dump-%6.6i"%i,"amuse",
                        append_to_file=False)    
    

def average(N,r,dat):
  n=len(r)
  a=numpy.argsort(r)
  i=0
  r_a=[]
  dat_a=[]
  while i < n:
    ra=r[a[i:i+N]].sum()/min(n-i,N)
    da=dat[a[i:i+N]].sum()/min(n-i,N)
    r_a.append(ra)
    dat_a.append(da)
    i=i+N
  return numpy.array(r_a),numpy.array(dat_a)


# below follow some plotting functions

def aplot(i, tag, xyc,xlim=None,ylim=None):
  pyplot.figure(figsize=(6,6))
  for x,y,c in xyc:
    xa,ya=average(100,x,y)
    pyplot.semilogy(xa,ya,c)  
  if xlim is not None:
    pyplot.xlim(xlim)
  if ylim is not None:
    pyplot.ylim(ylim)
  pyplot.savefig(tag+'-%6.6i.png'%i)

def xion_from_u(u):
  xion=xtrans*(u-uinit)/(utrans-uinit)
  a=numpy.where(u>utrans)[0]
  xion[a]=xtrans+(1-xtrans)*(u[a]-utrans)/(uion-utrans)
  return xion

def plots(i):

  g=read_set_from_file('dump-%6.6i'%i,'amuse')

  r=((g.x**2+g.y**2+g.z**2)**0.5).value_in(units.kpc)
  v=((g.vx**2+g.vy**2+g.vz**2)**0.5).value_in(units.kms)
  cs=(g.u**0.5).value_in(units.kms)
  xion=xion_from_u(g.u).value_in(units.none)
  rho=3./4/numpy.pi*8*g.mass/g.radius**3
  dens=(rho).value_in(units.amu/units.cm**3)
  pres=(g.u*rho).value_in(units.g/units.cm/units.s**2)
  mach=v/cs

  pyplot.figure(figsize=(6,6))
  pyplot.semilogy(r/15,xion,'r .')
  pyplot.semilogy(r/15,1-xion,'g .')
  pyplot.xlim((0.,1.))
  pyplot.ylim((1.e-6,1.))
  pyplot.savefig('xion-part-%6.6i.png'%i)

  aplot(i,'xion',((r/15,xion,'r'),(r/15,1-xion,'g')),
          xlim=(0.,1.),ylim=(1.e-6,1.))
  aplot(i,'pres',((r/15,pres,'r'),),
          xlim=(0.,1.),ylim=(1.e-17,1.e-14))
  aplot(i,'rho',((r/15,dens,'r'),),
          xlim=(0.,1.),ylim=(0.0001,0.01))
  aplot(i,'mach',((r/15,mach,'r'),),
          xlim=(0.,1.),ylim=(1.e-5,10.))

# main example

if __name__=="__main__":
  
  N=10000
  Ns=1
  L=15. | units.kpc
  dt=1. | units.Myr  
  tend=100.| units.Myr

  sph,rad=iliev_test_5(N,Ns,L,dt/2.)

  radhydro_evolve(sph,rad,tend,dt)
  if IS_PLOT_AVAILABLE:
    plots(50)
    plots(100)

