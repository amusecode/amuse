"""
   Nbody integration of N particles in N-body units from t=0 to
   t_end=1 Myr.  The initial cluster is a King (1966) model with
   dimension-less depth of the potential of W0=7. The initial
   distribution of stars is in virial equilibrium.  At this moment a
   4th order Hermite integrator is used for the integration.  Stellar
   masses are selected randomly from a Salpeter initial mass function
   between a minimum mass of Mmin=0.1MSun and Mmax=100MSun.  In order
   to assure proper scaling to astrophysical units, we have to define
   the cluster radius in physical units, in this case, we opted for a
   virial radius of 1pc.
   cluster in orbit in static potential
"""
import math
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse.units.optparse import OptionParser

class MilkyWay_galaxy(object):

    def __init__(self, Mb=1.40592e10| units.MSun, Md=8.5608e10| units.MSun, Mh=1.07068e11 | units.MSun  ):
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
        pot_disk= -constants.G*self.Md/(R**2+ (a2+ (z**2+ b2**2)**0.5 )**2 )**0.5
        #halo
        a3= 12.0 |units.kpc
        cut_off=100 |units.kpc
        d1= r/a3
        c=1+ (cut_off/a3)**1.02
        pot_halo= -constants.G*(self.Mh/a3)*d1**1.02/(1+ d1**1.02) -(constants.G*self.Mh/(1.02*a3))*(-1.02/c +numpy.log(c) + 1.02/(1+d1**1.02)- numpy.log(1.0 +d1**1.02) )
        return 2*(pot_bulge+pot_disk+ pot_halo) #I have to multiply by 2 because is a rigid potential

       
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
        az= force_bulge*z + force_disk*d*z/(z**2 + b2**2)**0.5  + force_halo*z/r 

        return ax,ay,az

    def vel_circ(self, r ):
        z=0 | units.kpc 
        b1= 0.3873 |units.kpc
        a2= 5.31 |units.kpc
        b2= 0.25 |units.kpc
        a3= 12.0 |units.kpc

        rdphi_b= constants.G*self.Mb*r**2/(r**2+b1**2)**1.5
        rdphi_d= constants.G*self.Md*r**2/(r**2+ (a2+(z**2+b2**2)**0.5)**2 )**1.5
        rdphi_h= constants.G*self.Mh*(r/a3)**0.02*r/(a3**2*(1+(r/a3)**1.02))

        vel_circb=  rdphi_b
        vel_circd= rdphi_d
        vel_circh= rdphi_h

        return (vel_circb+ vel_circd+ vel_circh)**0.5 

def main(N, W0, t_end, n_steps, filename, Mtot, Rvir, rgc, vgc):
    numpy.random.seed(111)

#    pos = (-1390, 9340, 25.3) |units.parsec
#    vel = (-207, -48.2, -6.72) | units.kms
    pSun = (8500, 0, 0) |units.parsec
    vSun_now = (-10.1, 235.5, 7.5) | units.kms
    vSun_back = (10.1, -235.5, -7.5) | units.kms

    pos = (6.091e+18,  2.508e+20, -9.487e+17) | units.m
    vel = (-2.460e+05, 9.160e+03, 7.090e+03) | units.ms

#    pos = [rgc.value_in(units.parsec),0,0] | units.parsec
#    vel = [0,vgc.value_in(units.kms),0] | units.kms

    converter=nbody_system.nbody_to_si(Mtot,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
#    bodies = new_plummer_model(N, W0,convert_nbody=converter)
#    bodies = Particles(1)
#    bodies.mass = 1 | units.MSun
#    bodies.position = (0,0,0)|units.kpc
#    bodies.velocity = (0,0,0)|units.kms
    eps2 = 0.25*(float(N))**(-0.666667) * Rvir**2
    bodies.scale_to_standard(convert_nbody=converter, 
                             smoothing_length_squared = eps2)
#    bodies.scale_to_standard(convert_nbody=converter)
    bodies.position += pos
    bodies.velocity += vel
    bodies.radius = 0 |  units.parsec

    cluster_gravity = ph4(converter)
    cluster_gravity.parameters.timestep_parameter = 0.01
    cluster_gravity.parameters.epsilon_squared = eps2
    cluster_gravity.particles.add_particles(bodies)
    channel_from_gravity_to_framework = cluster_gravity.particles.new_channel_to(bodies)
    
    write_set_to_file(bodies.savepoint(0.0|units.Myr), filename, 'hdf5')
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (MilkyWay_galaxy(),) )

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)
    t_orb = 2*numpy.pi*pos.length()/vel.length()
    gravity.timestep = min(dt/4., t_orb/32.)
    while time < t_end:
        time += dt
        gravity.evolve_model(time)

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_gravity_to_framework.copy()
        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

    print gravity.particles
    gravity.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 48,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-f", dest="filename", default = "proto_solar_cluster.hdf5",
                      help="output filename [%default]")
    result.add_option("-N", dest="N", type="int",default = 100,
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

