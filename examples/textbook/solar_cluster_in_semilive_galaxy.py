import math
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse.units.optparse import OptionParser
from amuse.units import quantities

from amuse.community.galaxia.interface import BarAndSpirals3D
from amuse.ext.composition_methods import *
from matplotlib import pyplot
from prepare_figure import figure_frame, get_distinct

class drift_without_gravity(object):

    def __init__(self, particles, time= 0 |units.Myr):
        self.particles=particles
        self.model_time= time
        self.softening_lengths_squared = (100|units.parsec)**2
        self.gravity_constant = constants.G
    def evolve_model(self, t_end):
        dt= t_end- self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time= t_end
    def add_particles(self, p):
        self.particles.add_particles(p)
    @property
    def potential_energy(self):
        return quantities.zero
    def get_potential_at_point(self,radius,x,y,z):
        positions = self.particles.position
        result = quantities.AdaptingVectorQuantity()
        for i in range(len(x)):
            dx = x[i] - positions.x
            dy = y[i] - positions.y
            dz = z[i] - positions.z
            dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
            dr = (dr_squared + self.softening_lengths_squared).sqrt()
            energy_of_this_particle = (self.particles.mass / dr).sum()
            result.append(-self.gravity_constant * energy_of_this_particle)
        return result
    def get_gravity_at_point(self,radius,x,y,z):
        positions = self.particles.position
        m1 = self.particles.mass
        result_ax = quantities.AdaptingVectorQuantity()
        result_ay = quantities.AdaptingVectorQuantity()
        result_az = quantities.AdaptingVectorQuantity()
        for i in range(len(x)):
            dx = x[i] - positions.x
            dy = y[i] - positions.y
            dz = z[i] - positions.z
            dr_squared = ((dx * dx) + (dy * dy) + (dz * dz) +
                self.softening_lengths_squared + radius[i]**2)

            ax = -self.gravity_constant * (m1*dx/dr_squared**1.5).sum()
            ay = -self.gravity_constant * (m1*dy/dr_squared**1.5).sum()
            az = -self.gravity_constant * (m1*dz/dr_squared**1.5).sum()

            result_ax.append(ax)
            result_ay.append(ay)
            result_az.append(az)
        return result_ax, result_ay, result_az
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()
    def stop(self):
        return

class MilkyWay_galaxy(object):

    def __init__(self, Mb=1.40592e10| units.MSun,
                 Md=8.5608e10| units.MSun, Mh=1.07068e11 | units.MSun  ):
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
        pot_disk = -constants.G*self.Md/(R**2+(a2+(z**2+b2**2)**0.5)**2)**0.5
        #halo
        a3= 12.0 |units.kpc
        cut_off=100 |units.kpc
        d1= r/a3
        c=1+ (cut_off/a3)**1.02
        pot_halo = -constants.G*(self.Mh/a3)*d1**1.02/(1 + d1**1.02) \
                     - (constants.G*self.Mh/(1.02*a3))\
                         *(-1.02/c + numpy.log(c) \
                           + 1.02/(1+d1**1.02) \
                           - numpy.log(1.0 +d1**1.02))
        return 2*(pot_bulge + pot_disk + pot_halo) # multiply by 2 for
    						   # rigid potential

       
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
       
        ax= force_bulge*x + force_disk*x + force_halo*x/r
        ay= force_bulge*y + force_disk*y + force_halo*y/r
        az= force_bulge*z + force_disk*d*z/(z**2 + b2**2)**0.5 + force_halo*z/r 

        return ax,ay,az

    def vel_circ(self, r ):
        z=0 | units.kpc 
        b1= 0.3873 |units.kpc
        a2= 5.31 |units.kpc
        b2= 0.25 |units.kpc
        a3= 12.0 |units.kpc

        rdphi_b= constants.G*self.Mb*r**2 \
                   / (r**2+b1**2)**1.5
        rdphi_d= constants.G*self.Md*r**2 \
                   / (r**2+ (a2+(z**2+b2**2)**0.5)**2 )**1.5
        rdphi_h= constants.G*self.Mh*(r/a3)**0.02*r \
                   / (a3**2*(1+(r/a3)**1.02))

        vel_circb = rdphi_b
        vel_circd = rdphi_d
        vel_circh = rdphi_h

        return (vel_circb+ vel_circd+ vel_circh)**0.5 

    def stop(self):
        return

def plot(stars, GMCs):
    figure = figure_frame("X [kpc]", "Y [kpc]", xsize=8, ysize=8)
    colors = get_distinct(2)

    print numpy.mean(stars.mass.value_in(units.MSun))
    print numpy.mean(GMCs.mass.value_in(units.MSun))
    size = stars.mass/(0.1 |units.MSun)
    pyplot.scatter(stars.x.value_in(units.kpc), stars.y.value_in(units.kpc),
                   s=size, c=colors[0])
    size = numpy.sqrt(GMCs.mass/(100 |units.MSun))
    pyplot.scatter(GMCs.x.value_in(units.kpc), GMCs.y.value_in(units.kpc),
                   s=size, alpha =  0.5, lw=0, c=colors[1])

    pyplot.scatter([0], [0], marker="+", s=300, c='r')
    pyplot.scatter([-8.4], [0], marker="o", s=100, c='g')    
    pyplot.axis("equal")
    pyplot.xlim(-10, 10)
    pyplot.ylim(-10, 10)
    pyplot.show()

def evolve_cluster_in_potential(gravity, t_end, dt,
                                channels_to_framework, sun=None, GMCs=None,
                                filename=False):
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    if filename:
        write_set_to_file(sun, filename, "hdf5", timestamp = time,
                          append_to_file=False)
        write_set_to_file(GMCs, filename, "hdf5", timestamp = time)

    x = []
    y = []
    while time < t_end:
        time += dt

        gravity.evolve_model(time)
        for ch in channels_to_framework:
            ch.copy()
        if filename:
            write_set_to_file(sun, filename, "hdf5", timestamp = time)
            write_set_to_file(GMCs, filename, "hdf5", timestamp = time)
        
        x.append(gravity.particles[0].x.value_in(units.kpc))
        y.append(gravity.particles[0].y.value_in(units.kpc))

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

    return x, y

def integrate_single_particle_in_potential(sun, t_end, dt):
    MWG = MilkyWay_galaxy()    
    cluster_gravity = drift_without_gravity(sun)
    channel_from_gravity_to_framework \
        = cluster_gravity.particles.new_channel_to(sun)
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (MWG,) )
    t_orb = 2*numpy.pi*sun.position.length()/sun.velocity.length()
    gravity.timestep = min(dt, 10|units.Myr)

    x, y = evolve_cluster_in_potential(gravity, t_end, dt,
                                       [channel_from_gravity_to_framework])
    gravity.stop()
    return x, y

def integrate_cluster_and_GMCs_in_potential(sun, GMCs, t_end, dt, filename):
    MWG = MilkyWay_galaxy()    
    GMC_gravity = drift_without_gravity(GMCs)
    channels = []
    channels.append(GMC_gravity.particles.new_channel_to(GMCs))
    
    converter = nbody_system.nbody_to_si(sun.mass.sum(),
                                         sun[0].position.length())
    cluster_gravity = BHTree(converter)
    cluster_gravity.particles.add_particles(sun)
    channels.append(cluster_gravity.particles.new_channel_to(sun))

    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (MWG, GMC_gravity) )
    gravity.add_system(GMC_gravity, (MWG,) )
    gravity.timestep = min(dt, 10|units.Myr)
    
    t_orb = 2*numpy.pi*sun.position.length()/sun.velocity.length()
    gravity.timestep = min(dt, 10|units.Myr)
    
    x, y = evolve_cluster_in_potential(gravity, t_end, dt,
                                       channels, sun, GMCs, filename)
    gravity.stop()
    return x, y

def initialize_sun_in_milky_way():
    sun = Particles(1)
    sun.mass= 1 | units.MSun
    sun.radius= 1 |units.RSun
    sun.position= [-8400.0, 0.0, 17.0] | units.parsec
    MWG = MilkyWay_galaxy()    
    vc = MWG.vel_circ(sun.position.length())
    sun.velocity= [11.352, (12.24+vc.value_in(units.kms)), 7.41] | units.kms
    current_sun = sun.copy()
    sun.velocity *= -1
    print "current:", sun
    return sun

def make_giant_molecular_clouds(Ngmc):
    N_thick_disk = int(0.5*Ngmc)
    N_thin_disk = int(0.5*Ngmc)
    converter=nbody_system.nbody_to_si(1.e+8|units.MSun, 1.0|units.kpc)
    from amuse.ext.protodisk import ProtoPlanetaryDisk
    Rin = 3.5 | units.kpc
    Rout = 7.5 | units.kpc
    masses = new_powerlaw_mass_distribution(N_thick_disk, alpha=-1.6,
                                            mass_min=1.0e+3|units.MSun,
                                            mass_max=1.0e+8|units.MSun)
    MGMCs = masses.sum()
    MWG = MilkyWay_galaxy()    
    v_inner = MWG.vel_circ(Rout)
    MGalaxy = v_inner**2*Rout/constants.G
    print "Masses:", MGMCs.in_(units.MSun), MGalaxy.in_(units.MSun), \
          MGMCs/MGalaxy
    GMCs = ProtoPlanetaryDisk(len(masses), convert_nbody=converter,
                              Rmin=Rin.value_in(units.kpc), 
                              Rmax=Rout.value_in(units.kpc),
                              q_out=30.0, discfraction=MGMCs/MGalaxy).result

    #second population of GMCs
    masses = new_powerlaw_mass_distribution(len(GMCs), alpha=-1.6,
                                            mass_min=1.e+3|units.MSun,
                                            mass_max=1.0e+8|units.MSun)    
    GMCs.mass = masses
    MGMCs = masses.sum()
    thin_disk_GMCs = ProtoPlanetaryDisk(N_thin_disk, convert_nbody=converter,
                              Rmin=Rin.value_in(units.kpc), 
                              Rmax=2*Rout.value_in(units.kpc),
                              q_out=10.0, discfraction=MGMCs/MGalaxy).result
    thin_disk_GMCs.masses = masses
    GMCs.add_particles(thin_disk_GMCs)
    GMCs.velocity *= -1
    GMCs.mass = new_powerlaw_mass_distribution(len(GMCs), alpha=-1.6,
                                                 mass_min=1.e+3|units.MSun,
                                                 mass_max=1.0e+8|units.MSun)
    print "v=", v_inner.in_(units.kms)
    print "GMC mass=", GMCs.mass.sum().in_(units.MSun)
    for gi in range(len(GMCs)):
        r = GMCs[gi].position.length()
        vc = MWG.vel_circ(r)
        GMCs[gi].velocity = GMCs[gi].velocity * (vc/GMCs[gi].velocity.length())

    return GMCs

def make_new_cluster(Ncl, Rvir, W0, sun):
    masses = new_salpeter_mass_distribution(Ncl, 0.1|units.MSun,
                                            10.0|units.MSun)
    converter = nbody_system.nbody_to_si(masses.sum(), Rvir)
    cluster = new_king_model(len(masses), W0=3, convert_nbody=converter)
    cluster.mass = masses
    eps2 = 0.25*len(masses)**(-2./3.) * Rvir**2
    cluster.scale_to_standard(convert_nbody=converter,
                              smoothing_length_squared = eps2)
    cluster.position += sun.position
    cluster.velocity += sun.velocity
    cluster.radius = 0 |  units.AU
    return cluster

def main(Ngmc, Ncl, W0, t_end, n_steps, filename, Rvir):
    numpy.random.seed(111)
    dt = t_end/float(n_steps)

    sun = initialize_sun_in_milky_way()

    print "Find birth location of the Sun."
    x, y = integrate_single_particle_in_potential(sun, t_end, dt)
    sun.velocity *= -1
    print "Birth location of the Sun:", sun

    GMCs = make_giant_molecular_clouds(Ngmc)
    cluster = make_new_cluster(Ncl, Rvir, W0, sun)

    integrate_cluster_and_GMCs_in_potential(cluster, GMCs, t_end, dt,
                                            filename)
    
    plot(cluster, GMCs)

def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 200,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-f", dest="filename",
                      default = "proto_solar_cluster.hdf5",
                      help="output filename [%default]")
    result.add_option("--Ncl", dest="Ncl", type="int",default = 1000,
                      help="number of stars [%default]")
    result.add_option("--Ngmc", dest="Ngmc", type="int",default = 1000,
                      help="number of GMCs [%default]")
    result.add_option("-R", unit= units.parsec,
                      dest="Rvir", type="float",default = 100 | units.parsec,
                      help="cluser virial radius [%default]")
    result.add_option("-t", unit= units.Gyr,
                      dest="t_end", type="float", default = 4.6 | units.Gyr,
                      help="end time of the simulation [%default]")
    result.add_option("-W", 
                      dest="W0", type="float", default = 7.0,
                      help="Dimensionless King potential depth (W0) [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

