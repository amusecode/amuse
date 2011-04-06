"""
  example of molecular cloud evolution with explictly 
  split SPH and grav evolution
  
  
"""  

import numpy
  
from matplotlib import pyplot 

from amuse.support.units import nbody_system
from amuse.support.units import units
    
from amuse.community.fi.interface import Fi

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.ext.derived_grav_systems import copycat
from amuse.ext.bridge import bridge

nc=2
nr=2
nplot=nc*nr
  
def make_map(sph,fig,t,i,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.parsec(x)
    y=units.parsec(y)
    z=units.parsec(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)

    rho=rho.reshape((N+1,N+1))

    s=fig.add_subplot(nc,nr,i)
    s.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
        extent=[-L/2,L/2,-L/2,L/2],vmin=2,vmax=6)
  
  
def run_mc(N=5000,Mcloud=10000. | units.MSun,Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

    dt=0.005 | units.Myr
    tend=0.12 | units.Myr

    parts=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube).result
    parts.h_smooth=0. | units.parsec

    sph=Fi(conv)

    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=False
    sph.parameters.gamma=1
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=dt/2  
    sph.parameters.verbosity=0 | units.none

    sph.gas_particles.add_particles(parts)

    grav=copycat(Fi, sph, conv)

    sys=bridge(verbose=False)
    sys.add_system(sph,(grav,),False)

    fig=pyplot.figure(figsize=(12,12))

    i=1
    make_map(sph,fig,0.*tend,i,N=200,L=3)
    while i<nplot:
      i=i+1
      sys.evolve_model((i-1)*tend/(nplot-1),timestep=dt) 
      make_map(sph,fig,(i-1)*tend/(nplot-1),i,N=200,L=3)

    pyplot.show()
#    pyplot.savefig('test.png')
  
if __name__ in ["__main__","__plot__"]:
    run_mc(10000,Mcloud=10000. | units.MSun,Rcloud=1. | units.parsec)
  
