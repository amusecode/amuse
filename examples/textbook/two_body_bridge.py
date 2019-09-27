import numpy 
from amuse.lab import *
from amuse.couple import bridge
from amuse.ext.solarsystem import new_solar_system

def main():
    filename = "SunAndEarth.hdf"
    ss = new_solar_system()
    star = ss[0]
    planet = ss[3]
    converter=nbody_system.nbody_to_si(star.mass,planet.position.length())
###BOOKLISTSTART###
    star_gravity = ph4(converter)
    star_gravity.particles.add_particle(star)

    planet_gravity = ph4(converter)
    planet_gravity.particles.add_particle(planet)

    channel_from_star_to_framework = star_gravity.particles.new_channel_to(ss)
    channel_from_planet_to_framework = planet_gravity.particles.new_channel_to(ss)

    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(star_gravity, (planet_gravity,) )
    gravity.add_system(planet_gravity, (star_gravity,) )
###BOOKLISTSTOP###

    write_set_to_file(ss, filename, 'hdf5', append_to_file=False)

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    gravity.timestep = 1|units.day
    time = zero
    dt = 10|units.day
    t_end = 100 | units.yr
    while time < t_end:
        time += dt
        gravity.evolve_model(time)

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_star_to_framework.copy()
        channel_from_planet_to_framework.copy()
        write_set_to_file(ss, filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print("T=", time, end=' ') 
        print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot) 
        Etot_prev = Etot
    gravity.stop()

if __name__ in ('__main__', '__plot__'):
    main()

