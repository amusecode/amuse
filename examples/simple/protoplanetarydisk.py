"""
  simple example of setting up a protoplanetary disk around a sun-like star

"""


import numpy

from matplotlib import pyplot

from amuse.community.fi.interface import Fi

from amuse.units import units
from amuse.units import nbody_system

from amuse.ext.protodisk import ProtoPlanetaryDisk

from amuse.datamodel import Particles
def make_map(sph,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.AU(x)
    y=units.AU(y)
    z=units.AU(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))

    return numpy.transpose(rho)
                    
if __name__ in ("__main__","__plot__"):

    N=20000
    tend=1. | units.yr
    Mstar=1. | units.MSun
        
    convert=nbody_system.nbody_to_si(Mstar, 1. | units.AU)
    proto=ProtoPlanetaryDisk(N,convert_nbody=convert,densitypower=1.5,Rmin=4,Rmax=20,q_out=1.)
    gas=proto.result
    gas.h_smooth=0.06 | units.AU 
             
    sun=Particles(1)
    sun.mass=Mstar
    sun.radius=2. | units.AU
    sun.x=0.|units.AU
    sun.y=0.|units.AU
    sun.z=0.|units.AU
    sun.vx=0.|units.kms
    sun.vy=0.|units.kms
    sun.vz=0.|units.kms
 
    sph=Fi(convert)
 
    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=True
    sph.parameters.gamma=1.
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=0.125 | units.yr  

    sph.gas_particles.add_particles(gas)
    sph.particles.add_particles(sun)
        
    sph.evolve_model(tend)    
            
    L=50
    print "1"
    rho=make_map(sph,N=200,L=L)
    sph.stop()
    pyplot.figure(figsize=(8,8))
    pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
        extent=[-L/2,L/2,-L/2,L/2],vmin=10,vmax=15)    
    pyplot.title(tend)
    pyplot.xlabel('AU')
    pyplot.savefig('test.png')
#    pyplot.show()
         
