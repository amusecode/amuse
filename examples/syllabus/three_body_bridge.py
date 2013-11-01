import numpy 
from amuse.lab import *
from amuse.couple import bridge

def new_system_of_sun_and_earth():
    stars = Particles(3)
    sun = stars[0]
    sun.mass = units.MSun(1.0)
    sun.position = units.m(numpy.array((0.0,0.0,0.0)))
    sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
    sun.radius = units.RSun(1.0)

    earth = stars[1]
    earth.mass = units.kg(5.9736e24)
    earth.radius = units.km(6371) 
    earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
    earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))

    moon = stars[2]
    moon.mass = units.kg(7.3477e22 )
    moon.radius = units.km(1737.10) 
    moon.position = units.km(numpy.array((149.5e6 + 384399.0 ,0.0,0.0)))
    moon.velocity = ([0.0,1.022,0] | units.km/units.s) + earth.velocity    
    return stars

def main():
    filename = "SunAndEarthAndMoon.hdf"
    ss = new_system_of_sun_and_earth()
    star = ss[0]
    planet = ss[1]
    moon = ss[2]
    converter=nbody_system.nbody_to_si(star.mass,planet.position.length())
    star_gravity = ph4(converter)
    star_gravity.particles.add_particle(star)

    planet_gravity = ph4(converter)
    planet_gravity.particles.add_particle(planet)

    moon_gravity = ph4(converter)
    moon_gravity.particles.add_particle(moon)

    channel_from_star_to_framework = star_gravity.particles.new_channel_to(ss)
    channel_from_planet_to_framework = planet_gravity.particles.new_channel_to(ss)
    channel_from_moon_to_framework = moon_gravity.particles.new_channel_to(ss)

    write_set_to_file(ss, filename, 'hdf5')
    
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(star_gravity, (planet_gravity,moon_gravity) )
    gravity.add_system(planet_gravity, (star_gravity,moon_gravity) )
    gravity.add_system(moon_gravity, (star_gravity,planet_gravity) )

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    gravity.timestep = 1|units.day
    time = zero
    dt = 3|units.day
    t_end = 1 | units.yr
    while time < t_end:
        time += dt
        gravity.evolve_model(time)

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_star_to_framework.copy()
        channel_from_planet_to_framework.copy()
        channel_from_moon_to_framework.copy()
        write_set_to_file(ss, filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot
    gravity.stop()

if __name__ in ('__main__', '__plot__'):
    main()

