from __future__ import print_function
import numpy
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot
from amuse.ic.kingmodel import new_king_model

"""
Implements a code simulating the galactic center. As the center itself does
not evolve we only need to define the 'get_gravity_at_point'
and 'get_potential_at_point'. Note that both functions get arrays
of points.
"""
class GalacticCenterGravityCode(object):
    def __init__(self,R, M, alpha):
        self.radius=R
        self.mass=M
        self.alpha=alpha

    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        m=self.mass*(r/self.radius)**self.alpha  
        fr=constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az

    def circular_velocity(self,r):  
        m=self.mass*(r/self.radius)**self.alpha  
        vc=(constants.G*m/r)**0.5
        return vc

    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2)**0.5
        c=constant.G*self.mass/self.radius**self.alpha    
        phi=c/(alpha-1)*(r**(self.alpha-1)-R**(self.alpha-1))
        return phi    
        
def make_king_model_cluster(nbodycode, N, W0, Mcluster,
                            Rcluster, parameters = []):

    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    bodies=new_king_model(N,W0,convert_nbody=converter)

    code=nbodycode(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code

def plot_cluster(x, y):

    from prepare_figure import single_frame, get_distinct
    colors = get_distinct(1)     
    f = single_frame('X [pc]', 'Y [pc]')
    pyplot.xlim(-60, 60)
    pyplot.ylim(-60, 60)
        
    pyplot.scatter(x,y, c=colors[0], s=50, lw=0)
    pyplot.savefig("Arches")
#    pyplot.show()

def evolve_cluster_in_galaxy(N, W0, Rinit, tend, timestep, M, R):

    Rgal=1. | units.kpc
    Mgal=1.6e10 | units.MSun
    alpha=1.2
    galaxy_code=GalacticCenterGravityCode(Rgal, Mgal, alpha)

    cluster_code=make_king_model_cluster(BHTree,N,W0, M,R,
        parameters=[("epsilon_squared", (0.01 | units.parsec)**2)])
    
    stars=cluster_code.particles.copy()    
    stars.x += Rinit
    stars.vy = 0.8*galaxy_code.circular_velocity(Rinit)
    channel=stars.new_channel_to(cluster_code.particles)
    channel.copy_attributes(["x","y","z","vx","vy","vz"])

    system=bridge(verbose=False)
    system.add_system(cluster_code, (galaxy_code,))

    times=numpy.arange(0|units.Myr, tend, timestep)
    for i,t in enumerate(times):
        system.evolve_model(t,timestep=timestep)
          
    x=system.particles.x.value_in(units.parsec)
    y=system.particles.y.value_in(units.parsec)
    cluster_code.stop()
    return x, y

if __name__ == "__main__":
    N=1024
    W0=3
    Rinit=50. | units.parsec
    timestep=0.01 | units.Myr
    endtime = 2.5 | units.Myr
    Mcluster = 5.e4 | units.MSun
    Rcluster = 0.8 | units.parsec
    x, y = evolve_cluster_in_galaxy(N, W0, Rinit, endtime, timestep, Mcluster, Rcluster)
    plot_cluster(x, y)
