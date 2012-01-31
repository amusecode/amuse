"""
  example of molecular cloud evolution with explictly 
  split SPH and grav evolution

  Initial condition is a smooth spherical cloud with random velocities
  as in Bonnell et al. (2003)  
  
"""  

import numpy
  
from matplotlib import pyplot 

from amuse.units import nbody_system
from amuse.units import units

from amuse.community.fi.interface import Fi

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.ext.derived_grav_systems import copycat
from amuse.ext.bridge import bridge

def make_map(sph,N=100,L=1):

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

    return rho
    
def run_mc(N=5000,Mcloud=10000. | units.MSun,Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

    dt=0.01 | units.Myr
    tend=0.36 | units.Myr

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
    sph.parameters.verbosity = 0

    sph.gas_particles.add_particles(parts)

    grav=copycat(Fi, sph, conv)

    sys=bridge(verbose=False)
    sys.add_system(sph,(grav,),False)

    fig=pyplot.figure(figsize=(12,12))

    ncolumn=2
    nrow=2
    nplot=ncolumn*nrow

    i=0
    L=3
    while i<nplot:
        ttarget=i*tend/(nplot-1)
        sys.evolve_model(ttarget,timestep=dt) 
        rho=make_map(sph,N=200,L=L)
        subplot=fig.add_subplot(ncolumn,nrow,i+1)
        subplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
        extent=[-L/2,L/2,-L/2,L/2],vmin=1,vmax=5)
        subplot.set_title(ttarget.in_(units.Myr))
        i=i+1

    sph.stop()
    pyplot.show()
#    pyplot.savefig('test.png')
  
if __name__ in ("__main__","__plot__"):
    run_mc(10000,Mcloud=1000. | units.MSun,Rcloud=1. | units.parsec)
  
