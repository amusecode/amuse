import numpy
from amuse.lab import *
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.phigrape.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot
from amuse.ic.kingmodel import new_king_model

###BOOKLISTSTART1###
class drift_without_gravity(object):
    def __init__(self, convert_nbody, time= 0 |units.Myr):
        self.model_time = time
        self.convert_nbody = convert_nbody
        self.particles = Particles()
    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time = t_end
    @property
    def potential_energy(self):
        return quantities.zero
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass \
                   *self.particles.velocity.lengths()**2).sum()
    def stop(self):
        pass
###BOOKLISTSTOP1###

class MilkyWay_galaxy(object):

    def __init__(self, Mb=1.40592e10| units.MSun,
                 Md=8.5608e10| units.MSun,
                 Mh=1.07068e11 | units.MSun):
        self.Mb = Mb
        self.Md = Md
        self.Mh = Mh

    def get_potential_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2)**0.5
        R = (x**2+y**2)**0.5
        # bulge
        b1 = 0.3873 |units.kpc
        pot_bulge = -constants.G*self.Mb/(r**2+b1**2)**0.5 
        # disk
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        pot_disk = -constants.G*self.Md/(R**2+(a2+(z**2+b2**2)**0.5)**2)**0.5
        #halo
        a3 = 12.0 |units.kpc
        cut_off = 100 |units.kpc
        d1 =  r/a3
        c = 1+ (cut_off/a3)**1.02
        pot_halo = -constants.G*(self.Mh/a3)*d1**1.02/(1+ d1**1.02) \
                   - (constants.G*self.Mh/(1.02*a3))\
                     * (-1.02/c +numpy.log(c) + 1.02/(1+d1**1.02) \
                        - numpy.log(1.0+d1**1.02))
        return 2*(pot_bulge+pot_disk+ pot_halo) # multiply by 2 for
    						# a rigid potential

    def get_gravity_at_point(self, eps, x,y,z): 
        r = (x**2+y**2+z**2)**0.5
        R = (x**2+y**2)**0.5
        #bulge
        b1 = 0.3873 |units.kpc
        force_bulge = -constants.G*self.Mb/(r**2+b1**2)**1.5 
        #disk
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        d = a2+ (z**2+ b2**2)**0.5
        force_disk =-constants.G*self.Md/(R**2+ d**2 )**1.5
        #halo
        a3 = 12.0 |units.kpc
        d1 = r/a3
        force_halo = -constants.G*self.Mh*d1**0.02/(a3**2*(1+d1**1.02))
       
        ax = force_bulge*x + force_disk*x + force_halo*x/r
        ay = force_bulge*y + force_disk*y + force_halo*y/r
        az = force_bulge*z + force_disk*d*z/(z**2 + b2**2)**0.5 \
                           + force_halo*z/r 

        return ax,ay,az

    def circular_velocity(self,r):  
        z = 0 | units.kpc 
        b1 = 0.3873 |units.kpc
        a2 = 5.31 |units.kpc
        b2 = 0.25 |units.kpc
        a3 = 12.0 |units.kpc

        rdphi_b = constants.G*self.Mb*r**2/(r**2+b1**2)**1.5
        rdphi_d = constants.G*self.Md*r**2/(r**2+(a2+(z**2+b2**2)**0.5)**2)**1.5
        rdphi_h = constants.G*self.Mh*(r/a3)**0.02*r/(a3**2*(1+(r/a3)**1.02))

        vel_circb = rdphi_b
        vel_circd = rdphi_d
        vel_circh = rdphi_h

        return (vel_circb + vel_circd + vel_circh)**0.5 
        
def make_king_model_cluster(nbodycode, N, W0, Mcluster,
                            Rcluster, parameters = []):

    converter = nbody_system.nbody_to_si(Mcluster,Rcluster)
    bodies = new_king_model(N,W0,convert_nbody=converter)

    code = nbodycode(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code

def plot_cluster(x, y):

    from prepare_figure import single_frame, get_distinct
    colors = get_distinct(1)     
    f = single_frame('X [kpc]', 'Y [kpc]')
    pyplot.xlim(-10, 10)
    pyplot.ylim(-10, 10)
        
    pyplot.scatter(x,y, c=colors[0], s=50, lw=0)
    pyplot.show()

def evolve_cluster_in_galaxy(N, W0, Rinit, tend, timestep, M, R):

    galaxy_code = MilkyWay_galaxy()

    converter = nbody_system.nbody_to_si(M, R)
    cluster_code = drift_without_gravity(convert_nbody=converter)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    cluster_code.particles.add_particles(bodies)
    
    stars = cluster_code.particles.copy()    
    stars.x += Rinit
    stars.vy = 0.8*galaxy_code.circular_velocity(Rinit)
    channel = stars.new_channel_to(cluster_code.particles)
    channel.copy_attributes(["x","y","z","vx","vy","vz"])

    system = bridge(verbose=False)
    system.add_system(cluster_code, (galaxy_code,))

    times = quantities.arange(0|units.Myr, tend, 100*timestep)
    for i,t in enumerate(times):
        print("Time=", t.in_(units.Myr))
        system.evolve_model(t, timestep=timestep)
          
    x = system.particles.x.value_in(units.kpc)
    y = system.particles.y.value_in(units.kpc)
    cluster_code.stop()
    return x, y

if __name__ == "__main__":
    N = 1024
    W0 = 3
    Rinit = 8.5 | units.kpc
    timestep = 0.1 | units.Myr
    endtime = 4.6 | units.Gyr
    Mcluster = 5.e4 | units.MSun
    Rcluster = 1.0 | units.parsec
    x, y = evolve_cluster_in_galaxy(N, W0, Rinit, endtime, timestep,
                                    Mcluster, Rcluster)
    plot_cluster(x, y)
