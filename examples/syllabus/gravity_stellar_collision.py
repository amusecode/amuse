"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
from __future__ import print_function
from amuse.lab import *
from amuse.io import store
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


    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_se_to_framework.copy_attributes(["mass","radius", "age"])
    
    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5')
    E_init = gravity.kinetic_energy + gravity.potential_energy

    
    Nenc = 0
    #dE_enc = dE_dyn = dE_stellar = zero
    dE_coll = zero
    time = zero
    while time < t_end:
        time += dt

        bodies.radius *= 1.e+5
        channel_from_framework_to_gd.copy_attributes(["mass", "radius"])
        E_dyn = gravity.kinetic_energy  + gravity.potential_energy 
        gravity.evolve_model(time)
        dE_dyn = E_dyn - (gravity.kinetic_energy  + gravity.potential_energy)

        if stopping_condition.is_set():
            E_coll = gravity.kinetic_energy + gravity.potential_energy
            print("At time=", gravity.model_time.in_(units.Myr), "number of encounters=", len(stopping_condition.particles(0)))
            for ci in range(len(stopping_condition.particles(0))): 
                particles_in_encounter = Particles(particles=[stopping_condition.particles(0)[ci], stopping_condition.particles(1)[ci]])
                particles_in_encounter = particles_in_encounter.get_intersecting_subset_in(bodies)

                merge_two_stars(bodies, particles_in_encounter)
                bodies.synchronize_to(gravity.particles)
                bodies.synchronize_to(stellar.particles)
                Nenc+=1
                print("Resolve encounter Number:", Nenc)
                gravity.evolve_model(time)
            dE_coll = E_coll - (gravity.kinetic_energy + gravity.potential_energy)

        E_stellar = gravity.kinetic_energy + gravity.potential_energy 
        stellar.evolve_model(time)
        dE_stellar = E_stellar - (gravity.kinetic_energy + gravity.potential_energy)

        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["age", "mass", "radius"])

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')
        print_diagnostics(time, bodies.mass.sum(), E_dyn, dE_dyn, dE_coll, dE_stellar)

    gravity.stop()
    stellar.stop()

def print_diagnostics(time, Mtot, Etot, dE_dyn, dE_coll, dE_stellar):
        print("T=", time, end=' ') 
        print("M=", Mtot, end=' ') 
        print("E= ", Etot, end=' ') 
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
                      help="maximal stellar mass [%default.value_in(units.MSun))]")
    result.add_option("-m", unit=units.MSun,
                      dest="Mmin", type="float",default = 0.1,
                      help="minimal stellar mass [%default.value_in(units.MSun)]")
    result.add_option("-R", unit=units.parsec,
                      dest="Rvir", type="float",default = 1.0,
                      help="cluser virial radius [%default.value_in(units.parsec)]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 10.0,
                      help="end time of the simulation [%default.value_in(units.Myr]")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [%default]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

