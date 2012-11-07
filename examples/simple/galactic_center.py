"""
Evolves a cluster in the potention of the galactic center
  
Uses the bridge integrator to couple different codes.
In this example a cluster is evolved  circling the galactic center, represented by a static potential.
"""

import numpy

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system

from amuse.ext.bridge import bridge
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2

from matplotlib import pyplot

from amuse.ic.kingmodel import new_king_model

class GalacticCenterGravityCode(object):
    """
    Implements a code simulating the galactic center. As the center itself does
    not evolve we only need to define the 'get_gravity_at_point'
    and 'get_potential_at_point'. Note that both functions get arrays
    of points.
    """
    def __init__(self,R=1000.| units.parsec, M=1.6e10 | units.MSun, alpha=1.2):
        self.R=R
        self.M=M
        self.alpha=alpha

    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        m=self.M*(r/self.R)**self.alpha  
        fr=constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az

    def get_potential_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        c=constant.G*self.M/self.R**self.alpha    
        phi=c/(alpha-1)*(r**(self.alpha-1)-R**(self.alpha-1))
        return phi    

    def vcirc(self,r):  
        m=self.M*(r/self.R)**self.alpha  
        vc=(constants.G*m/r)**0.5
        return vc
        
def king_model_cluster(interface,N=1024,W0=3, Mcluster=4.e4 | units.MSun,
                                 Rcluster= .7 | units.parsec,parameters=[]):
    """
    helper function to setup an nbody king model cluster (returns code with particles)
    """      
    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)

    parts=new_king_model(N,W0,convert_nbody=converter)
    parts.radius=0.0| units.parsec

    nb=interface(converter)
    for name,value in parameters:
      setattr(nb.parameters, name, value)

    nb.particles.add_particles(parts)

    return nb

def shift_sys(system,dx,dy,dz,dvx,dvy,dvz):
    """
    helper function to shift system
    """
    parts=system.particles.copy()
    parts.x=parts.x+dx
    parts.y=parts.y+dy
    parts.z=parts.z+dz
    parts.vx=parts.vx+dvx
    parts.vy=parts.vy+dvy
    parts.vz=parts.vz+dvz
    channel=parts.new_channel_to(system.particles)
    channel.copy_attributes(["x","y","z","vx","vy","vz"])
    #      parts.copy_values_of_state_attributes_to(system.particles)
    
if __name__ in ('__main__', '__plot__'):

# parameter setup:
    N=1024
    W0=3
    Rinit=50. | units.parsec
    timestep=0.01 | units.Myr
    Mcluster=4.e4 | units.MSun
    Rcluster=0.7 | units.parsec

# make cluster and Galactic center
    cluster=king_model_cluster(Fi,N,W0, Mcluster,Rcluster, parameters=[
                   ("epsilon_squared", (0.01 | units.parsec)**2), 
                   ("periodic_box_size",200 | units.parsec),
                   ("timestep",timestep/4)] )
    center=GalacticCenterGravityCode()
    
# shift the cluster to an orbit around GC
# (note the systems share the same coordinate frame, although units may differ)    
    vcirc=center.vcirc(Rinit)
    shift_sys(cluster,Rinit,0| units.parsec,0.|units.parsec,
                  0| units.kms,0.8*vcirc,0| units.kms)

# setup bridge; cluster is evolved under influence of GC
    sys=bridge(verbose=False)
    sys.add_system(cluster, (center,), False)   


# evolve and make plots
    times=units.Myr([0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4])
    f=pyplot.figure(figsize=(8,16))

    for i,t in enumerate(times):
        sys.evolve_model(t,timestep=timestep)

        x=sys.particles.x.value_in(units.parsec)
        y=sys.particles.y.value_in(units.parsec)

        subplot=f.add_subplot(4,2,i+1)
        subplot.plot(x,y,'r .')
        subplot.plot([0.],[0.],'b +')
        subplot.set_xlim(-60,60)
        subplot.set_ylim(-60,60)
        subplot.set_title(t)
        if i==7:
            subplot.set_xlabel('parsec')
            
    cluster.stop()
    pyplot.show()
#    pyplot.savefig('test.eps')
    
