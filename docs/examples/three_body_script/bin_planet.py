import numpy
import time
from matplotlib import pyplot

from amuse.units import nbody_system,units,constants
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.io import write_set_to_file,read_set_from_file

from amuse.community.huayno.interface import Huayno


from amuse.datamodel import Particles

def binary(m1=1.|units.MSun,m2=1.| units.MSun,r1=None,r2=None,ecc=0,P=1| units.yr):
  mu=constants.G*(m1+m2)
  a=(P/(2*numpy.pi)*mu**0.5)**(2./3.)

  f1=m2/(m1+m2)
  f2=m1/(m1+m2)
  rmax=a*(1+ecc)

  r0=rmax
    
  h=(a*mu*(1-ecc**2))**0.5
  v0=h/r0

  bin=Particles(2)

  bin[0].mass=m1
  bin[0].x=r0*f1
  bin[0].vy=v0*f1
  bin[1].mass=m2
  bin[1].x=-r0*f2
  bin[1].vy=-v0*f2

  bin.y=0*r0
  bin.z=0.*r0
  bin.vx=0*v0
  bin.vz=0.*v0
  if r1 is None:
    bin[0].radius=(1.|units.RSun)*(m1/(1.|units.MSun))**(1./3.)
  else:
    bin[0].radius=r1
  if r2 is None:
    bin[1].radius=(1.|units.RSun)*(m2/(1.|units.MSun))**(1./3.)
  else:
    bin[1].radius=r2

  return bin

def binary_with_planet(m1=1.|units.MSun, m2=1.| units.MSun, m_planet=1|units.MJupiter,
                       r1=None, r2=None, r_planet=None, ecc_binary=0, P_binary=20 | units.day,
                         ecc_planet=0., P_planet=1.| units.yr, pangle_planet=0., a_planet=None):

  parts=binary(m1,m2,r1,r2,ecc_binary,P_binary)


  mu=constants.G*(m1+m2+m_planet)
  if a_planet is None:
    if P_planet is None:
      print "provide a_planet or P_planet"
      raise Exception  
    a_planet=(P_planet/(2*numpy.pi)*mu**0.5)**(2./3.)
      
  rmax=a_planet*(1+ecc_planet)
  r0=rmax

  print a_planet
  print a_planet.in_(units.AU),r0.in_(units.AU)
  
  h=(a_planet*mu*(1-ecc_planet**2))**0.5
  v0=h/r0

  planet=Particles(1)
  planet.mass=m_planet
  if r_planet is None:
    r_planet=(1.|units.RJupiter)*(m_planet/(1.|units.MJupiter))**(1./3.)
  planet.radius=r_planet
  planet.x=numpy.cos(pangle_planet)*r0
  planet.y=numpy.sin(pangle_planet)*r0
  planet.z=0.*r0
  planet.vx=-numpy.sin(pangle_planet)*v0
  planet.vy=numpy.cos(pangle_planet)*v0
  planet.vz=0.*v0

  parts.add_particles(planet)

  parts.move_to_center()

  return parts
  

def elements(starmass,x,y,z,vx,vy,vz,G=constants.G):
    mu=G*starmass
    r=(x**2+y**2+z**2)**0.5
    v2=(vx**2+vy**2+vz**2)
    
    e=v2/2-mu/r
    
    a=-mu/2/e
    
    hx=y*vz-z*vy
    hy=z*vx-x*vz
    hz=x*vy-y*vx

    rdotv=x*vx+y*vy+z*vz

    ex=v2*x/mu-rdotv*vx/mu-x/r
    ey=v2*y/mu-rdotv*vy/mu-y/r
    ez=v2*z/mu-rdotv*vz/mu-z/r

    pangle= numpy.arctan2(ex,ey) # assuming orbits in plane
    pangle=pangle

    h2=(hx**2+hy**2+hz**2)    
    
    eps=(1-h2/mu/a)**0.5
    
    return a,eps,pangle    
    
def binary_with_planet_run(m1=1.|units.MSun, m2=1.| units.MSun, m_planet=1|units.MJupiter,
                           r1=None, r2=None, r_planet=None, ecc_binary=0, P_binary=20 | units.day,
                           ecc_planet=0., P_planet=1.| units.yr, pangle_planet=0., a_planet=None,
                           tend=100. | units.yr,hostname=''):

    dEcrit=1.e-6

    three=binary_with_planet(m1=m1,m2=m2,m_planet=m_planet,r1=r1,r2=r2,r_planet=r_planet,
             ecc_binary=ecc_binary,P_binary=P_binary,
             ecc_planet=ecc_planet,a_planet=a_planet,pangle_planet=pangle_planet)
    
    convert=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
    code=Huayno(convert,hostname=hostname)
    
    code.parameters.inttype_parameter=code.inttypes.SHARED10
    code.parameters.timestep_parameter=.2
    code.parameters.timestep=100. | units.day
        
    dt=10000. | units.day
    
    code.particles.add_particles(three)

    E0=code.potential_energy+code.kinetic_energy
    a0,eps0,pangle0=elements( three.total_mass(),
          code.particles.x[2],
          code.particles.y[2],
          code.particles.z[2],
          code.particles.vx[2],
          code.particles.vy[2],
          code.particles.vz[2] )

    t=0. | units.day
    result="stable"
    while(t < tend-dt/2):
        t=t+dt
        code.evolve_model(t)
        E=code.potential_energy+code.kinetic_energy    
        dE=abs(((E-E0)/E0))
        a,eps,pangle=elements( three.total_mass(),
          code.particles.x[2],
          code.particles.y[2],
          code.particles.z[2],
          code.particles.vx[2],
          code.particles.vy[2],
          code.particles.vz[2] )
        if dE > dEcrit or a<0.5*a0 or a>2.*a0:
          result="unstable"
          if dE > dEcrit:
            result="failed"
          break
    code.stop()
    return result,t,dE,a.in_(units.AU),eps,pangle
    
def test_run():        
    three=binary_with_planet(
      m1=0.6897 | units.MSun,m2=0.20255 | units.MSun,m_planet=0.333 | units.MJupiter,
      r1=0.6489 | units.RSun,r2=0.22623 | units.RSun,r_planet=0.754 | units.RJupiter,
      ecc_binary=0.15944,P_binary=41.08| units.day,ecc_planet=0.00685,a_planet=.7048 | units.AU,
      pangle_planet=0.) 

    convert=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
    code=Huayno(convert)
    
    code.parameters.inttype_parameter=code.inttypes.SHARED4
    code.parameters.timestep_parameter=0.1
    
#    tend=100. | units.yr
    tend=100. | units.day
    snapfreq=1
    
    dt=10. | units.day
#    dt=convert.to_si( 1. | nbody_system.time).in_(units.day)
    
    code.particles.add_particles(three)

    x = AdaptingVectorQuantity()
    y = AdaptingVectorQuantity()
    z = AdaptingVectorQuantity()
    vx = AdaptingVectorQuantity()
    vy = AdaptingVectorQuantity()
    vz = AdaptingVectorQuantity()
    x.append(code.particles.x)
    y.append(code.particles.y)
    z.append(code.particles.z)
    vx.append(code.particles.vx)
    vy.append(code.particles.vy)
    vz.append(code.particles.vz)
    ts=[0.]
    E0=code.potential_energy+code.kinetic_energy    
    dE=[1.e-14]
    t=0. | units.day
    i=0
    while(t < tend-dt/2):
        i+=1
        t=t+dt
        if i%snapfreq==0:
          print t
          ts.append(t.value_in(units.day))
          code.evolve_model(t)
          x.append(code.particles.x)
          y.append(code.particles.y)
          z.append(code.particles.z)
          vx.append(code.particles.vx)
          vy.append(code.particles.vy)
          vz.append(code.particles.vz)
          E=code.potential_energy+code.kinetic_energy    
          dE.append(abs(((E-E0)/E0)))
    code.stop()

    a,eps,pangle=elements(three.total_mass(),
    x[:,2],
    y[:,2],
    z[:,2],
    vx[:,2],
    vy[:,2],
    vz[:,2])

    x=x.value_in(units.AU)
    y=y.value_in(units.AU)
    
    a=a.value_in(units.AU)
    eps=eps
    
    print a[-1],eps[-1],pangle[-1]

    f=pyplot.figure(figsize=(8,8))
    pyplot.plot(x[:,0],y[:,0],'r.')
    pyplot.plot(x[:,1],y[:,1],'g.')
    pyplot.plot(x[:,2],y[:,2],'b.')
    pyplot.xlim(-3,3)
    pyplot.ylim(-3,3)
    pyplot.xlabel('AU')
    pyplot.savefig('three_16b.eps')

    f=pyplot.figure(figsize=(8,8))
    pyplot.semilogy(ts,dE,'g.')
    pyplot.xlabel('time (day)')
    pyplot.ylabel('dE/E0')
    pyplot.savefig('three_16b_eerr.eps')

    f=pyplot.figure(figsize=(8,8))
    pyplot.plot(ts,a,'g.')
    pyplot.xlabel('time (day)')
    pyplot.ylabel('a (AU)')
    pyplot.savefig('three_16b_a.eps')

    f=pyplot.figure(figsize=(8,8))
    pyplot.plot(ts,eps,'g.')
    pyplot.xlabel('time (day)')
    pyplot.ylabel('eccentricity')
    pyplot.savefig('three_16b_ecc.eps')

    f=pyplot.figure(figsize=(8,8))
    pyplot.plot(ts,pangle,'g.')
    pyplot.xlabel('time (day)')
    pyplot.ylabel('long. of periapsis')
    pyplot.savefig('three_16b_pangle.eps')

if __name__=="__main__":
    result,t,dE,a,eps,pangle=binary_with_planet_run(
      m1=0.6897 | units.MSun,m2=0.20255 | units.MSun,m_planet=0.333 | units.MJupiter,
      r1=0.6489 | units.RSun,r2=0.22623 | units.RSun,r_planet=0.754 | units.RJupiter,
      ecc_binary=0.52,P_binary=41.08| units.day,ecc_planet=0.00685,a_planet=.87048 | units.AU,
      pangle_planet=0., tend=1000.| units.yr,hostname="gaasp") 
    print result,t,dE,a,eps,pangle

