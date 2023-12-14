import numpy
import struct

from .interface import Kepler

from amuse.units import nbody_system
from amuse.units import units,constants

from amuse.ic.plummer import new_plummer_model

from amuse.datamodel import Particle

#from matplotlib import pyplot

import time

from amuse.ext.orbital_elements import orbital_elements_from_binary,new_binary_from_orbital_elements

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

    h2=(hx**2+hy**2+hz**2)

    eps=(1-h2/mu/a)**0.5

    return a,eps

def test_kepler(N=10000, tend=1.| units.yr,method=0):

  numpy.random.seed(12345)

  conv=nbody_system.nbody_to_si(2.| units.MSun, 5.|units.AU)

  comets=new_plummer_model(N,conv)

  sun=Particle(mass=1.|units.MSun)

  sun.position=[0,0,0]|units.AU
  sun.velocity=[0,0,0]|units.kms

  comets.mass*=0.

  code=Kepler(conv,redirection="none")

  code.set_method(method)

  code.central_particle.add_particle(sun)
  code.orbiters.add_particles(comets)

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz)


#  print code.orbiters.x[0]
  print(orbital_elements_from_binary(code.particles[0:2],constants.G))

  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print(orbital_elements_from_binary(code.particles[0:2],constants.G))
#  print code.orbiters.x[0]



  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  dev=numpy.where(da > 0.00001)[0]
  print(len(dev))

  print(a0[dev].value_in(units.AU))
  print(eps0[dev])

#  pyplot.plot(a0[dev].value_in(units.AU),eps0[dev],"ro")
#  pyplot.plot(a[dev].value_in(units.AU),eps[dev],"g+")

  print("max da,deps:",da.max(), deps.max())

  print("time:",t2-t1)

#  pyplot.show()

  return t2-t1,da.max(),deps.max()

def test_kepler_almost_parabolic( tend=1,method=0):
  code=Kepler(redirection="none")

  code.set_method(method)

  mass1=1.| nbody_system.mass
  mass2=0| nbody_system.mass
  semimajor_axis=1.|nbody_system.length
  eccentricity=0.9999999
  p=2*numpy.pi*(semimajor_axis**3/nbody_system.G/mass1)**0.5
  tend=tend*p
  print(tend)
  parts=new_binary_from_orbital_elements(
          mass1,
          mass2,
          semimajor_axis,
          eccentricity = eccentricity,
          true_anomaly = 0.0102121
      )

  code.central_particle.add_particle(parts[0])
  code.orbiters.add_particle(parts[1])

  a0,eps0=elements(mass1,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  print(orbital_elements_from_binary(code.particles[0:2]))

  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print(orbital_elements_from_binary(code.particles[0:2]))

  print(code.orbiters.position)

  a,eps=elements(mass1,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print(da,deps)
  print("time:",t2-t1)

def test_kepler_parabolic( tend=1,method=0, sign=+1):
  code=Kepler(redirection="none")

  code.set_method(method)

  sun=Particle()
  sun.mass=1. | nbody_system.mass
  sun.x=0. | nbody_system.length
  sun.y=0. | nbody_system.length
  sun.z=0. | nbody_system.length
  sun.vx=0. | nbody_system.speed
  sun.vy=0. | nbody_system.speed
  sun.vz=0. | nbody_system.speed

  comet=Particle()
  comet.mass= 0 | nbody_system.mass
  comet.x=1. | nbody_system.length
  comet.y=0. | nbody_system.length
  comet.z=0. | nbody_system.length
  comet.vx=0. | nbody_system.speed
  comet.vy=(1.0 + sign * 1.0e-10)*(2*nbody_system.G*sun.mass/comet.x)**0.5
  comet.vz=0. | nbody_system.speed

  tend=tend | nbody_system.time
  print(tend)

  code.central_particle.add_particle(sun)
  code.orbiters.add_particle(comet)

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  print(orbital_elements_from_binary(code.particles[0:2]))

  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print(orbital_elements_from_binary(code.particles[0:2]))

  print(code.orbiters.position)

  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print(da,deps)
  print("time:",t2-t1)

def crash_test(method=1):
  code=Kepler(redirection="none")

  code.set_method(method)

  smu=1.224744871391589
  mu=smu**2
  r0=2.787802728537455
  rv0=-0.9899959571994231
  alpha=0.01380749549277993
  smudt=2.809925892593303
  v02=(mu*(2/r0-alpha))
  vx=rv0
  vy=(v02-vx**2)**0.5

  sun=Particle()
  sun.mass=mu | nbody_system.mass
  sun.x=0. | nbody_system.length
  sun.y=0. | nbody_system.length
  sun.z=0. | nbody_system.length
  sun.vx=0. | nbody_system.speed
  sun.vy=0. | nbody_system.speed
  sun.vz=0. | nbody_system.speed

  comet=Particle()
  comet.mass= 0 | nbody_system.mass
  comet.x=r0| nbody_system.length
  comet.y=0. | nbody_system.length
  comet.z=0. | nbody_system.length
  comet.vx=vx | nbody_system.speed
  comet.vy=vy | nbody_system.speed
  comet.vz=0. | nbody_system.speed

  tend=(smudt/smu) | nbody_system.time
  print(tend)

  code.central_particle.add_particle(sun)
  code.orbiters.add_particle(comet)

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  print(orbital_elements_from_binary(code.particles[0:2]))

  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print(orbital_elements_from_binary(code.particles[0:2]))

  print(code.orbiters.position)

  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print(da,deps)
  print("time:",t2-t1)


def crash_test2(method=1):
  code=Kepler(redirection="none")

  code.set_method(method)

  """
  mu=struct.unpack('!d','3ff7ffffffffffff'.decode('hex'))[0]
  dt=struct.unpack('!d','40025ab746b00001'.decode('hex'))[0]
  pos1=struct.unpack('!d','bfed36dc82998ed4'.decode('hex'))[0]
  pos2=struct.unpack('!d','40051297fc6e5256'.decode('hex'))[0]
  pos3=struct.unpack('!d','0000000000000000'.decode('hex'))[0]
  vel1=struct.unpack('!d','3fb09d8008ba33b9'.decode('hex'))[0]
  vel2=struct.unpack('!d','bff06788b551b81d'.decode('hex'))[0]
  vel3=struct.unpack('!d','0000000000000000'.decode('hex'))[0]
  """
  mu=float.fromhex("0x1.8p+0")
  dt=float.fromhex("0x1.25ab746bp+1")
  pos1=float.fromhex("-0x1.d36dc82998ed4p-1")
  pos2=float.fromhex("0x1.51297fc6e5256p+1")
  pos3=float.fromhex("0x0p+0")
  vel1=float.fromhex("0x1.09d8008ba33b9p-4")
  vel2=float.fromhex("-0x1.06788b551b81ep+0")
  vel3=float.fromhex("0x0p+0")

  sun=Particle()
  sun.mass=mu | nbody_system.mass
  sun.x=0. | nbody_system.length
  sun.y=0. | nbody_system.length
  sun.z=0. | nbody_system.length
  sun.vx=0. | nbody_system.speed
  sun.vy=0. | nbody_system.speed
  sun.vz=0. | nbody_system.speed

  comet=Particle()
  comet.mass= 0 | nbody_system.mass
  comet.x=pos1 | nbody_system.length
  comet.y=pos2 | nbody_system.length
  comet.z=pos3 | nbody_system.length
  comet.vx=vel1 | nbody_system.speed
  comet.vy=vel2 | nbody_system.speed
  comet.vz=vel3 | nbody_system.speed

  tend=dt | nbody_system.time
  print(tend,mu)

  code.central_particle.add_particle(sun)
  code.orbiters.add_particle(comet)

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  print(orbital_elements_from_binary(code.particles[0:2]))

  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print(orbital_elements_from_binary(code.particles[0:2]))

  print(code.orbiters.position)

  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print(da,deps)
  print("time:",t2-t1)

def test_softening(method=1):
  code=Kepler(redirection="none")

  code.set_method(method)

  dt=float.fromhex("0x1.67b39e372f04dp+4")
  mu=float.fromhex("0x1.fffffffffffdfp-3")
  e2=float.fromhex("0x1.0000000000003p+0")
  pos1=float.fromhex("0x1.1b76542265052p-1")
  pos2=float.fromhex("0x1.0c4dbda42097cp-6")
  pos3=float.fromhex("0x1.54fd66cd1e212p-3")
  vel1=float.fromhex("0x1.d6ef43d58ca7ep-2")
  vel2=float.fromhex("0x1.7a85379e59794p-2")
  vel3=float.fromhex("-0x1.5421044d1acffp-1")

  sun=Particle()
  sun.mass=mu | nbody_system.mass
  sun.x=0. | nbody_system.length
  sun.y=0. | nbody_system.length
  sun.z=0. | nbody_system.length
  sun.vx=0. | nbody_system.speed
  sun.vy=0. | nbody_system.speed
  sun.vz=0. | nbody_system.speed

  comet=Particle()
  comet.mass= 0 | nbody_system.mass
  comet.x=pos1 | nbody_system.length
  comet.y=pos2 | nbody_system.length
  comet.z=pos3 | nbody_system.length
  comet.vx=vel1 | nbody_system.speed
  comet.vy=vel2 | nbody_system.speed
  comet.vz=vel3 | nbody_system.speed

  tend=dt | nbody_system.time
  print(tend,mu)

  code.central_particle.add_particle(sun)
  code.orbiters.add_particle(comet)

  code.parameters.epsilon_squared = e2 | nbody_system.length**2

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  print(orbital_elements_from_binary(code.particles[0:2]))

  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print(orbital_elements_from_binary(code.particles[0:2]))

  print(code.orbiters.position)

  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print(da,deps)
  print("time:",t2-t1)

def t_linear(tend=1,N=100,method=0):
  code=Kepler(redirection="none")

  code.set_method(method)

  mass=1. | nbody_system.mass
  x=1. | nbody_system.length
  vx=0 | nbody_system.speed

  e=0.5*vx**2-nbody_system.G*mass/x

  semimajor_axis=-nbody_system.G*mass/2/e

  p=2*numpy.pi*(semimajor_axis**3/nbody_system.G/mass)**0.5

  print(semimajor_axis)
  print(p)

  tend=tend*p
  dt=p/N

  sun=Particle()
  sun.mass=mass
  sun.x=0. | nbody_system.length
  sun.y=0. | nbody_system.length
  sun.z=0. | nbody_system.length
  sun.vx=0. | nbody_system.speed
  sun.vy=0. | nbody_system.speed
  sun.vz=0. | nbody_system.speed

  comet=Particle()
  comet.mass= 0 | nbody_system.mass
  comet.x=x
  comet.y=0. | nbody_system.length
  comet.z=0. | nbody_system.length
  comet.vx=vx
  comet.vy=0. | nbody_system.speed
  comet.vz=0. | nbody_system.speed

  code.central_particle.add_particle(sun)
  code.orbiters.add_particle(comet)

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  print(orbital_elements_from_binary(code.particles[0:2]))

  #pyplot.ion()
  #f=pyplot.figure(figsize=(8,6))
  #pyplot.show()

  tnow=0*tend
  time=[]
  xs=[]
  while tnow<tend:
    tnow+=dt
    print(tnow,int(tnow/dt))
    code.evolve_model(tnow)
    #f.clf()
    time.append(tnow/tend)
    xs.append(code.orbiters.x[0].number)
  #pyplot.plot(time,xs,"r+")
  #pyplot.xlim(-0.1,1.1)
  #pyplot.ylim(-1.1,3.1)
  #pyplot.draw()

  print(orbital_elements_from_binary(code.particles[0:2]))

  print(code.orbiters.position)

  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)

  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print(da,deps)
  input()


if __name__=="__main__":

  for method in [1,0]:
      t_linear(N=100,method=method)
      print()

  print("-"*10)
  print()


  tend = 1.0

  for method in [1,0]:
      crash_test(method=method)
      print()

  print("-"*10)
  print()

  for method in [1,0]:
      crash_test2(method=method)
      print()

  print("-"*10)
  print()

  for method in [1,0]:
    test_kepler_parabolic(tend=tend,method=method, sign=+1)
    print()

  print("-"*10)
  print()

  for method in [1,0]:
    test_kepler_parabolic(tend=tend,method=method, sign=-1)
    print()

  print("-"*10)
  print()

  for method in [1,0]:
    test_kepler_almost_parabolic(tend=tend,method=method)
    print()

  print("-"*10)
  print()

  for method in [1,0]:
    test_kepler(N=10000,tend=tend | units.yr,method=method)
    print()

  print("-"*10)
  print()

  for method in [0,]:
      test_softening(method=method)
      print()

