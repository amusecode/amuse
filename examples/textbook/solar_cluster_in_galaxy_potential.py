import math
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse.units.optparse import OptionParser
from amuse.units import quantities

from amuse.community.galaxia.interface import BarAndSpirals3D
from amuse.ext.composition_methods import *
from matplotlib import pyplot

from prepare_figure import single_frame
from distinct_colours import get_distinct


class drift_without_gravity(object):

    def __init__(self, particles, time= 0 |units.Myr):
        self.particles=particles
        self.model_time= time
    def evolve_model(self, t_end):
        dt= t_end- self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time= t_end
    def add_particles(self, p):
        self.particles.add_particles(p)
    @property
    def potential_energy(self):
        return quantities.zero
    @property 
    def kinetic_energy(self):
        return \
            (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()
    def stop(self):
        return

class MilkyWay_galaxy(object):

    def __init__(self, Mb=1.40592e10| units.MSun,
                 Md=8.5608e10| units.MSun,
                 Mh=1.07068e11 | units.MSun  ):
        self.Mb= Mb
        self.Md= Md
        self.Mh= Mh

    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2)**0.5
        R= (x**2+y**2)**0.5
        # buldge
        b1= 0.3873 |units.kpc
        pot_bulge= -constants.G*self.Mb/(r**2+b1**2)**0.5 
        # disk
        a2= 5.31 |units.kpc
        b2= 0.25 |units.kpc
        pot_disk = \
            -constants.G*self.Md/(R**2 + (a2+ (z**2+ b2**2)**0.5 )**2 )**0.5
        #halo
        a3= 12.0 |units.kpc
        cut_off=100 |units.kpc
        d1= r/a3
        c=1+ (cut_off/a3)**1.02
        pot_halo= -constants.G*(self.Mh/a3)*d1**1.02/(1+ d1**1.02) \
                  - (constants.G*self.Mh/(1.02*a3))\
                      * (-1.02/c +numpy.log(c) + 1.02/(1+d1**1.02) \
                           - numpy.log(1.0 +d1**1.02) )
        return 2*(pot_bulge+pot_disk+ pot_halo) # multiply by 2 because it
    						# is a rigid potential

       
    def get_gravity_at_point(self, eps, x,y,z): 
        r= (x**2+y**2+z**2)**0.5
        R= (x**2+y**2)**0.5
        #bulge
        b1= 0.3873 |units.kpc
        force_bulge= -constants.G*self.Mb/(r**2+b1**2)**1.5 
        #disk
        a2= 5.31 |units.kpc
        b2= 0.25 |units.kpc
        d= a2+ (z**2+ b2**2)**0.5
        force_disk=-constants.G*self.Md/(R**2+ d**2 )**1.5
        #halo
        a3= 12.0 |units.kpc
        d1= r/a3
        force_halo= -constants.G*self.Mh*d1**0.02/(a3**2*(1+d1**1.02))
       
        ax= force_bulge*x + force_disk*x  + force_halo*x/r
        ay= force_bulge*y + force_disk*y  + force_halo*y/r
        az= force_bulge*z + force_disk*d*z/(z**2 + b2**2)**0.5 + force_halo*z/r 

        return ax,ay,az

    def vel_circ(self, r ):
        z=0 | units.kpc 
        b1= 0.3873 |units.kpc
        a2= 5.31 |units.kpc
        b2= 0.25 |units.kpc
        a3= 12.0 |units.kpc

        rdphi_b = constants.G*self.Mb*r**2/(r**2+b1**2)**1.5
        rdphi_d =constants.G*self.Md*r**2/(r**2+(a2+(z**2+b2**2)**0.5)**2 )**1.5
        rdphi_h = constants.G*self.Mh*(r/a3)**0.02*r/(a3**2*(1+(r/a3)**1.02))

        vel_circb =  rdphi_b
        vel_circd = rdphi_d
        vel_circh = rdphi_h

        return (vel_circb+ vel_circd+ vel_circh)**0.5 

    def stop(self):
        return
    
def evolve_cluster_in_potential(gravity, t_end, dt, channel_to_framework):
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    x = []
    y = []
    while time < t_end:
        time += dt

        gravity.evolve_model(time)
        channel_to_framework.copy()
        x.append(gravity.particles[0].x.value_in(units.kpc))
        y.append(gravity.particles[0].y.value_in(units.kpc))

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        """
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        """
        Etot_prev = Etot

    return x, y

def integrate_single_particle_in_potential(sun, t_end, dt, converter):

    # cluster_gravity = drift_without_gravity(sun)
    cluster_gravity = BHTree(converter)
    cluster_gravity.particles.add_particles(sun)
    channel_from_gravity_to_framework \
        = cluster_gravity.particles.new_channel_to(sun)
    
    MWG = MilkyWay_galaxy()
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (MWG,) )
    t_orb = 2*numpy.pi*sun.position.length()/sun.velocity.length()
    gravity.timestep = min(dt, 10|units.Myr)

    x, y = evolve_cluster_in_potential(gravity, t_end, dt,
                                       channel_from_gravity_to_framework)
    gravity.stop()
    return x, y

def main(N, W0, t_end, n_steps, filename, Mtot, Rvir, rgc, vgc):
    numpy.random.seed(111)
    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    dt = t_end/float(n_steps)

    sun = Particles(1)
    sun.mass= 1 | units.MSun
    sun.radius= 1 |units.RSun
    sun.position= [-8400.0, 0.0, 17.0] | units.parsec
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    fig = pyplot.figure(figsize=(8,8))	
    pyplot.xlim(-10, 10)
    pyplot.ylim(-10, 10)
    pyplot.axis('equal')
    pyplot.xlabel("X [kpc]")
    pyplot.ylabel("Y [kpc]")
    colors = get_distinct(6)

    MWG = MilkyWay_galaxy()    
    vc = MWG.vel_circ(sun.position.length())
    sun.velocity = [11.352, (12.24+vc.value_in(units.kms)), 7.41] | units.kms
    sun.velocity *= -1

    print("Current Sun:")
    print(sun)

    print("\nFinding birth location of the Sun...")
    x, y = integrate_single_particle_in_potential(sun, t_end, dt, converter)
    pyplot.plot(x, y, lw=4, alpha=0.2, c=colors[1])

    print("Initial Sun:")
    print(sun)
    sun.velocity *= -1

    cluster = new_king_model(N, W0=3, convert_nbody=converter)
    cluster.mass = new_salpeter_mass_distribution(len(cluster),
                                                  0.1|units.MSun,
                                                  10.0|units.MSun)
    eps2 = 0.25*(float(N))**(-0.666667) * Rvir**2
    cluster.scale_to_standard(convert_nbody=converter,
                              smoothing_length_squared = eps2)
    cluster.position += sun.position
    cluster.velocity += sun.velocity
    cluster.radius = 0 |  units.parsec
    
    pyplot.scatter(cluster.x.value_in(units.kpc),
                   cluster.y.value_in(units.kpc),
                   s=10, c=colors[3])
    print('\nTracking', N, 'siblings')
    x, y = integrate_single_particle_in_potential(cluster, t_end, dt,
                                                  converter)
    size = cluster.mass/(0.1 |units.MSun)
    pyplot.scatter(cluster.x.value_in(units.kpc),
                   cluster.y.value_in(units.kpc),
                   c=colors[0], alpha=1.0, lw=0, s=size)
    pyplot.scatter(sun.x.value_in(units.kpc),
                   sun.y.value_in(units.kpc),
                   marker='+', s=100, c=colors[2])

    pyplot.scatter([0], [0], marker="+", s=300, c='r')
    pyplot.scatter([-8.4], [0], marker="o", s=100, c='g')

    save_file = 'SolarClusterInPotential.png'
    pyplot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 2000,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-f", dest="filename",
                      default = "proto_solar_cluster.hdf5",
                      help="output filename [%default]")
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float",default = 100 | units.MSun,
                      help="cluster mass [%default]")
    result.add_option("-R", unit= units.parsec,
                      dest="Rvir", type="float",default = 100 | units.parsec,
                      help="cluser virial radius [%default]")
    result.add_option("-r", unit= units.parsec,
                      dest="rgc", type="float",default = 8500.0 | units.parsec,
                      help="distance to the galactic center [%default]")
    result.add_option("-v", unit= units.kms,
                      dest="vgc", type="float",default = 220.0 | units.kms,
            help="orbital velotiy around the galactic center [%default]")
    result.add_option("-t", unit= units.Gyr,
                      dest="t_end", type="float", default = 4.8 | units.Gyr,
                      help="end time of the simulation [%default]")
    result.add_option("-W", 
                      dest="W0", type="float", default = 7.0,
            help="Dimension-less depth of the King potential (W0) [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

