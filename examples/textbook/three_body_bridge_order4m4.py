import numpy 
from amuse.lab import *
from amuse.ext.composition_methods import *
from amuse.couple import bridge

from amuse.ext.solarsystem import new_solar_system

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
    filename = "SunAndEarthAndMoon_TBB.h5"
    ss = new_system_of_sun_and_earth()
    star = ss[0]
    planet = ss[1]
    moon = ss[2]

    print(moon)
    converter=nbody_system.nbody_to_si(star.mass, 1|units.AU)
    star_gravity = ph4(converter)
    star_gravity.particles.add_particle(star)

    planet_gravity = ph4(converter)
    planet_gravity.particles.add_particle(planet)

    moon_gravity = ph4(converter)
    moon_gravity.particles.add_particle(moon)

    channel_from_star_to_framework = star_gravity.particles.new_channel_to(ss)
    channel_from_planet_to_framework = planet_gravity.particles.new_channel_to(ss)
    channel_from_moon_to_framework = moon_gravity.particles.new_channel_to(ss)

    write_set_to_file(ss.savepoint(0.0|units.Myr), filename, 'hdf5', append_to_file=False, version='2')
    
    gravity = bridge.Bridge(use_threading=False, method=SPLIT_4TH_S_M4)
####gravity = bridge.Bridge(use_threading=False, method=SPLIT_6TH_SS_M13)
    gravity.add_system(star_gravity, (planet_gravity,moon_gravity) )
    gravity.add_system(planet_gravity, (star_gravity,moon_gravity) )
    gravity.add_system(moon_gravity, (star_gravity,planet_gravity) )

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init
    DDE_max = 0 

    gravity.timestep = 1|units.day
    time = zero
    dt = 30|units.day
    t_end = 100 | units.yr
    while time < t_end:
        time += dt
        gravity.evolve_model(time)

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_star_to_framework.copy()
        channel_from_planet_to_framework.copy()
        channel_from_moon_to_framework.copy()
        write_set_to_file(ss.savepoint(time), filename, 'hdf5', version='2')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        DDE = (Etot_prev-Etot)/Etot
        DDE_max = max(abs(DDE_max), abs(DDE))
        print("T=", time, end=' ') 
        print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", DDE, "dde_min=", DDE_max)
        Etot_prev = Etot
    gravity.stop()

if __name__ in ('__main__', '__plot__'):
    main()

