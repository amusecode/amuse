import numpy
from amuse.lab import *
from amuse.io import store
from amuse.community.seba.interface import SeBa

def merge_two_stars(bodies, particles_in_encounter):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle=Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.age = min(particles_in_encounter.age) \
                         * max(particles_in_encounter.mass)/new_particle.mass
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = 0 | units.RSun
    bodies.add_particles(new_particle)
    bodies.remove_particles(particles_in_encounter)

def resolve_collision(collision_detection, gravity, stellar, bodies):
    if collision_detection.is_set():
        E_coll = gravity.kinetic_energy + gravity.potential_energy
        print("At time=", gravity.model_time.in_(units.Myr), \
              "number of encounters=", len(collision_detection.particles(0)))
        Nenc = 0
        for ci in range(len(collision_detection.particles(0))): 
            particles_in_encounter \
                = Particles(particles=[collision_detection.particles(0)[ci],
                                       collision_detection.particles(1)[ci]])
            particles_in_encounter \
                = particles_in_encounter.get_intersecting_subset_in(bodies)

            merge_two_stars(bodies, particles_in_encounter)
            bodies.synchronize_to(gravity.particles)
            bodies.synchronize_to(stellar.particles)
            Nenc += 1
            print("Resolve encounter Number:", Nenc)
        dE_coll = E_coll - (gravity.kinetic_energy + gravity.potential_energy)
        print("dE_coll =", dE_coll, "N_enc=", Nenc)

###BOOKLISTSTART3###
def resolve_supernova(supernova_detection, bodies, time):
    if supernova_detection.is_set():
        print("At time=", time.in_(units.Myr), \
              len(supernova_detection.particles(0)), 'supernova(e) detected')

        Nsn = 0
        for ci in range(len(supernova_detection.particles(0))):
            print(supernova_detection.particles(0))
            particles_in_supernova \
                = Particles(particles=supernova_detection.particles(0))
            natal_kick_x = particles_in_supernova.natal_kick_x
            natal_kick_y = particles_in_supernova.natal_kick_y
            natal_kick_z = particles_in_supernova.natal_kick_z

            particles_in_supernova \
                = particles_in_supernova.get_intersecting_subset_in(bodies)
            particles_in_supernova.vx += natal_kick_x
            particles_in_supernova.vy += natal_kick_y
            particles_in_supernova.vz += natal_kick_z
            Nsn += 1

        print('Resolved', Nsn, 'supernova(e)')
###BOOKLISTSTOP3###
        
def main(N, W0, t_end, dt, filename, Rvir, Mmin, Mmax, z):

    numpy.random.seed(1)

###BOOKLISTSTART1###
    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(convert_nbody=converter)

    gravity = ph4(converter)
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)
    collision_detection = gravity.stopping_conditions.collision_detection
    collision_detection.enable()

    stellar = SeBa()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)
    stellar.evolve_model(0|units.Myr)

    supernova_detection = stellar.stopping_conditions.supernova_detection
    supernova_detection.enable()

    channel_from_se = stellar.particles.new_channel_to(bodies)
    channel_from_gd = gravity.particles.new_channel_to(bodies)
    channel_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_se.copy_attributes(["mass","radius", "age",
                                     "temperature", "luminosity"])
###BOOKLISTSTOP1###
    
    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5')
    E_init = gravity.kinetic_energy + gravity.potential_energy
    
    Nsn = 0
    dE_coll = zero
    time = zero

###BOOKLISTSTART2###
    while time < t_end:

        dt_min = min(dt, stellar.particles.time_step.min())
        print("Time steps:", dt.in_(units.Myr), dt_min.in_(units.Myr))
        time += dt_min

        bodies.radius *= 1.e+5
        channel_to_gd.copy_attributes(["mass", "radius", "vx", "vy", "vz"])
        E_dyn = gravity.kinetic_energy  + gravity.potential_energy 
        gravity.evolve_model(time)
        dE_dyn = E_dyn - (gravity.kinetic_energy  + gravity.potential_energy)

        resolve_collision(collision_detection, gravity, stellar, bodies)
        channel_from_gd.copy()

        time = gravity.model_time
        
        E_stellar = gravity.kinetic_energy + gravity.potential_energy 
        stellar.evolve_model(time)
        dE_stellar = E_stellar - (gravity.kinetic_energy
                                   + gravity.potential_energy)

        resolve_supernova(supernova_detection, bodies, time)
        channel_from_se.copy_attributes(["mass", "radius", "age",
                                         "temperature", "luminosity"])

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')
        print_diagnostics(time, bodies.mass.sum(), E_dyn,
                          dE_dyn, dE_coll, dE_stellar)
###BOOKLISTSTOP2###

    gravity.stop()
    stellar.stop()

def print_diagnostics(time, Mtot, Etot, dE_dyn, dE_coll, dE_stellar):
    print("T=", time.in_(units.Myr), end=' ') 
    print("M=", Mtot.in_(units.MSun), end=' ') 
    print("E= ", Etot, end=' ') 
    print("dE(dyn)=", dE_dyn/Etot, end=' ') 
    print("dE(coll)=", dE_coll/Etot, end=' ') 
    print("dE(se)=", dE_stellar/Etot)
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity_stellar.hdf5",
                      help="output filename [%default]")
    result.add_option("-N", dest="N", type="int",default = 10,
                      help="number of stars [%default]")
    result.add_option("--dt", unit=units.Myr,
                      dest="dt", type="float",default = 1|units.Myr,
                      help="output timesteps [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [%default.value_in(units.MSun))]")
    result.add_option("-m", unit=units.MSun,
                      dest="Mmin", type="float",default = 10,
                      help="minimal stellar mass [%default.value_in(units.MSun)]")
    result.add_option("-R", unit=units.parsec,
                      dest="Rvir", type="float",default = 0.1,
                      help="cluser virial radius [%default.value_in(units.parsec)]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 50.0,
                      help="end time of the simulation [%default.value_in(units.Myr]")
    result.add_option("-W", dest="W0", type="float", default = 12.0,
                      help="Dimension-less depth of the King potential (W0) [%default]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

