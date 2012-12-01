"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from amuse.community.seba.interface import SeBa

def merge_two_stars(bodies, particles_in_encounter):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle=Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.age = min(particles_in_encounter.age) * max(particles_in_encounter.mass)/new_particle.mass
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = 0 | units.RSun
    bodies.add_particles(new_particle)
    bodies.remove_particles(particles_in_encounter)
    
def main(N=10, W0=7.0, t_end=10, nsteps=10, filename="gravity_stellar.hdf5", Rvir=1, Mmin=0.1, Mmax= 100, z=0.02):
    t_end = t_end | nbody_system.time
    dt = t_end/float(nsteps)
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses

#    stellar = SeBa()
    stellar = SSE()
    stellar.parameters.metallicity = z
    bodies.radius = 1 * Rvir 
    stellar.particles.add_particle(bodies)
    stellar.commit_particles()

    bodies.radius = stellar.particles.radius

    gravity = ph4(converter)
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])

    bodies.scale_to_standard(convert_nbody=converter)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()
    
    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5')
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    t_end = converter.to_si(t_end)
    time = 0.0 | t_end.unit
    dt = converter.to_si(dt)
    Nenc = 0
    dEk_enc = zero    
    dEp_enc = zero
    while time < t_end:
        time += dt

        gravity.evolve_model(time)
        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        if stopping_condition.is_set():
            Ek_enc = gravity.kinetic_energy 
            Ep_enc = gravity.potential_energy
            print "At time=", time, "number of encounters=", len(stopping_condition.particles(0))
            for ci in range(len(stopping_condition.particles(0))): 
                particles_in_encounter = Particles(particles=[stopping_condition.particles(0)[ci], stopping_condition.particles(1)[ci]])
                particles_in_encounter = particles_in_encounter.get_intersecting_subset_in(bodies)

                merge_two_stars(bodies, particles_in_encounter)
                bodies.synchronize_to(gravity.particles)
                bodies.synchronize_to(stellar.particles)
                Nenc+=1
                print "Resolve encounter Number:", Nenc
            dEk_enc += Ek_enc - gravity.kinetic_energy 
            dEp_enc += Ep_enc - gravity.potential_energy

        stellar.evolve_model(time)

        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["age", "mass","radius","luminosity"])
        bodies.radius = 1.e+5*bodies.radius
        channel_from_framework_to_gd.copy_attributes(["mass", "radius"])
        bodies.radius = bodies.radius/1.e+5

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_prev-Etot
        dE_se = Etot_prev_se-Etot
        Mtot = bodies.mass.sum()
        print "T=", time, 
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, 
        print "(dE[SE]=", dE_se/Etot, ")"
        print "dE(enc)=", dEk_enc, dEp_enc
        Etot_init -= dE
        Etot_prev = Etot

    print bodies
    print len(bodies)
    print bodies.age

    gravity.stop()
    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity_stellar.hdf5",
                      help="output filename [gravity_stellar.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [100]")
    result.add_option("-n", dest="nsteps", type="int",default = 10,
                      help="number of output steps [10]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [100] MSun")
    result.add_option("-m", dest="Mmin", type="float",default = 0.1,
                      help="minimal stellar mass [0.1] MSun")
    result.add_option("-R", dest="Rvir", type="float",default = 1.0,
                      help="cluser virial radius [1] in parsec")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] Myr")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

