"""
   N-body integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metallicity z.
"""
from amuse.lab import *
from amuse.io import store
from amuse.community.seba.interface import SeBa

###BOOKLISTSTART1###
def merge_two_stars(bodies, particles_in_encounter):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle=Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.age = min(particles_in_encounter.age) \
                         * max(particles_in_encounter.mass)/new_particle.mass
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = particles_in_encounter.radius.sum()
    bodies.add_particles(new_particle)
    print("Two stars (M=",particles_in_encounter.mass, \
          ") collided with d=", com_pos.length())
    bodies.remove_particles(particles_in_encounter)
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def resolve_collision(collision_detection, gravity, stellar, bodies):
    if collision_detection.is_set():
        E_coll = gravity.kinetic_energy + gravity.potential_energy
        print("Collision at time=", gravity.model_time.in_(units.Myr))
        for ci in range(len(collision_detection.particles(0))): 
            particles_in_encounter \
                = Particles(particles=[collision_detection.particles(0)[ci],
                                       collision_detection.particles(1)[ci]])
            particles_in_encounter \
                = particles_in_encounter.get_intersecting_subset_in(bodies)
            d = (particles_in_encounter[0].position \
                   - particles_in_encounter[1].position).length()
            if particles_in_encounter.collision_radius.sum() > d:
                merge_two_stars(bodies, particles_in_encounter)
                bodies.synchronize_to(gravity.particles)
                bodies.synchronize_to(stellar.particles)
            else:
                print("Failed to resolve encounter: stars too small.")
            dE_coll = E_coll - (gravity.kinetic_energy \
                                  + gravity.potential_energy)
        print("Energy error in the collision: dE =", dE_coll) 
###BOOKLISTSTOP2###

def main(N, W0, t_end, dt, filename, Rvir, Mmin, Mmax, z):

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(convert_nbody=converter)

    gravity = ph4(converter)
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    stellar = SeBa()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)

    channel_from_stellar = stellar.particles.new_channel_to(bodies,
        attributes=["mass", "radius", "age"],
        target_names=["mass", "radius", "age"])
    channel_from_gravity = gravity.particles.new_channel_to(bodies,
        attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "radius"],
        target_names=["x", "y", "z", "vx", "vy", "vz",
                      "mass", "collision_radius"])
    channel_to_gravity = bodies.new_channel_to(gravity.particles,
                                    attributes=["mass", "collision_radius"],
                                    target_names=["mass", "radius"])
    channel_from_stellar.copy()
    
    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5')
    E_init = gravity.kinetic_energy + gravity.potential_energy

    Nenc = 0
    dE_coll = zero
    time = zero

###BOOKLISTSTART3###
    while time < t_end:
        time += dt

        E_stellar = gravity.kinetic_energy + gravity.potential_energy 
        stellar.evolve_model(time)
        dE_stellar = E_stellar - (gravity.kinetic_energy \
                                   + gravity.potential_energy)

        channel_from_stellar.copy()
        bodies.collision_radius = 1.e+5 * bodies.radius
        channel_to_gravity.copy()
        
        E_dyn = gravity.kinetic_energy  + gravity.potential_energy 
        gravity.evolve_model(time)
        dE_dyn = E_dyn - (gravity.kinetic_energy  + gravity.potential_energy)

        channel_from_gravity.copy()
        resolve_collision(stopping_condition, gravity, stellar, bodies)

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')
        print_diagnostics(time, bodies.mass.sum(),
                          E_dyn, dE_dyn, dE_coll, dE_stellar)
###BOOKLISTSTOP3###

    gravity.stop()
    stellar.stop()

def print_diagnostics(time, Mtot, Etot, dE_dyn, dE_coll, dE_stellar):
        print("Time=", time, end=' ') 
        print("Mtot=", Mtot, end=' ') 
        print("Etot= ", Etot, end=' ') 
        print("dE(dyn)=", dE_dyn/Etot, end=' ') 
        print("dE(coll)=", dE_coll/Etot, end=' ') 
        print("dE(se)=", dE_stellar/Etot)
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity_stellar.hdf5",
                      help="output filename [%default]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [%default]")
    result.add_option("--dt", unit=units.Myr,
                      dest="dt", type="float",default = 1|units.Myr,
                      help="output timesteps [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mmax", type="float",default = 100,
                help="maximum stellar mass [%default.value_in(units.MSun))]")
    result.add_option("-m", unit=units.MSun,
                      dest="Mmin", type="float",default = 0.1,
                help="minimum stellar mass [%default.value_in(units.MSun)]")
    result.add_option("-R", unit=units.parsec,
                      dest="Rvir", type="float",default = 1.0,
                help="cluser virial radius [%default.value_in(units.parsec)]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 10.0,
                help="end time of the simulation [%default.value_in(units.Myr]")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                help="Dimensionless depth King potential depth (W0) [%default]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    set_printing_strategy("custom", 
                          preferred_units=[units.MSun, units.RSun, units.Myr], 
                          precision=4, prefix = "", 
                          separator=" [", suffix="]")
    main(**o.__dict__)

