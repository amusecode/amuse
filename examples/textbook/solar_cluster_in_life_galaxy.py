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
from optparse import OptionParser
from amuse.units.optparse import OptionParser
from amuse.couple import bridge
from amuse.community.halogen.interface import HalogenInterface, Halogen
import amuse.community.halogen

def initialize_galaxy_model(N, converter):
    #default_options = dict(redirection = "none")
    instance = Halogen(converter)
    instance.initialize_code()
    instance.parameters.alpha = 2.0
    instance.parameters.beta  = 5.0
    instance.parameters.gamma = 0.0
    instance.parameters.number_of_particles = N
    instance.parameters.random_seed = 1
    instance.commit_parameters()
    instance.generate_particles()
    galaxy = instance.particles.copy()
        
    instance.cleanup_code()
    instance.stop()
    return galaxy

def main(N=1000, n=100, W0=7.0, t_end=10|units.Myr, dt_diag=1|units.Myr, filename="nbody.hdf5", Mtot=100|units.MSun, Rvir=1|units.parsec):

    pos = [-1.72566460731, -11.0998897304, -0.0790965379539] | units.kpc
    vel = [174.619387484, -36.7939713616, 0.85020621627] | units.kms

    masses = new_salpeter_mass_distribution(n, 0.1|units.MSun, 10.0|units.MSun)
    converter=nbody_system.nbody_to_si(masses.sum(),Rvir)
    cluster = new_king_model(n, W0,convert_nbody=converter)
    cluster.mass = masses
    cluster.scale_to_standard(convert_nbody=converter)
    cluster.position += pos
    cluster.velocity += vel
    cluster.radius = 0 |  units.parsec

    cluster_gravity = ph4(converter)
    cluster_gravity.parameters.timestep_parameter = 0.1
    cluster_gravity.parameters.epsilon_squared = (100|units.AU)**2
    cluster_gravity.particles.add_particles(cluster)
    channel_from_cluster_to_framework = cluster_gravity.particles.new_channel_to(cluster)

    stellar = SeBa()
    stellar.particles.add_particles(cluster)
    channel_from_stellar_to_cluster_gravity = stellar.particles.new_channel_to(cluster_gravity.particles)
    
    converter=nbody_system.nbody_to_si(2.067352e+11|units.MSun,1|units.kpc)
    galaxy = initialize_galaxy_model(N=N, converter=converter)
    print len(galaxy)
    galaxy_gravity = BHTree(converter)
    galaxy_gravity.parameters.timestep_parameter = 0.1
    galaxy_gravity.parameters.epsilon_squared = (100|units.parsec)**2
    galaxy_gravity.particles.add_particles(galaxy)
    channel_from_galaxy_to_framework = galaxy_gravity.particles.new_channel_to(galaxy)
    
    moving_bodies = ParticlesSuperset([cluster, galaxy])

    time = 0.0 | t_end.unit
    write_set_to_file(galaxy.savepoint(time), filename, 'hdf5', append_to_file=False)
    write_set_to_file(cluster.savepoint(time), filename, 'hdf5')
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (galaxy_gravity,) )
    gravity.add_system(galaxy_gravity, (cluster_gravity,) )

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    dt = dt_diag
    t_orb = 2*numpy.pi*pos.length()/vel.length()
    gravity.timestep = t_orb/32.
    while time < t_end:
        time += dt

        stellar.evolve_model(time)
        channel_from_stellar_to_cluster_gravity.copy()

        gravity.evolve_model(time)

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_cluster_to_framework.copy()
        channel_from_galaxy_to_framework.copy()
        write_set_to_file(galaxy.savepoint(time), filename, 'hdf5')
        write_set_to_file(cluster.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

    gravity.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-d", unit=units.Myr, 
                      dest="dt_diag", type="float", default = 0.5|units.Gyr,
                      help="diagnostics output time steps [%default]")
    result.add_option("-f", dest="filename", default = "proto_solar_cluster_life.hdf5",
                      help="output filename [%default]")
    result.add_option("-N", dest="N", type="int",default = 10000,
                      help="number of stars in the Galaxy [%default]")
    result.add_option("-n", dest="n", type="int",default = 1000,
                      help="number of stars in the cluster [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float",default = 10000|units.MSun,
                      help="cluster mass [100] MSun")
    result.add_option("-r", unit=units.parsec,
                      dest="Rvir", type="float",default = 100.0|units.parsec,
                      help="cluser virial radius [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 4.6|units.Gyr,
                      help="end time of the simulation [1] %unit")
    result.add_option("-W", dest="W0", type="float", default = 3.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

