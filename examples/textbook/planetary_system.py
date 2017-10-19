import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from optparse import OptionParser
from amuse.couple import bridge
import make_planets_oligarch
from amuse.ext.orbital_elements import orbital_elements_from_binary

def calculate_orbital_elements(star, planet):

    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, outc_out, lan_out, aop_out = orbital_elements_from_binary(p, G=constants.G)
    return a, e

from amuse.ext.orbital_elements import new_binary_from_orbital_elements
def make_planetesimal_disk(Nplanetesimals):

    amin = 6.4|units.AU
    amax = 30|units.AU
    emax = 0.1
    imax = 1
    a = amin + (amax-amin)*numpy.random.random_sample(Nplanetesimals)
    e = emax*numpy.random.random_sample(Nplanetesimals)
    i = imax*numpy.random.random_sample(Nplanetesimals)
    ta = 0 #numpy.acos(np.random.uniform(0,2*numpy.pi,Nplanetesimals))
    loan = 0
    aof = 0
    mp = 0.1 | units.MEarth
    planetesimals = Particles(Nplanetesimals)
    for i, pi in enumerate(planetesimals):
        b = new_binary_from_orbital_elements(Mstar, mp, a[i], e[i], ta, inc[i],
                                             loan, aof, G=constant.G)
        pi.mass = mp
        pi.position = b[1].position
        pi.velocity = b[1].velocity

    return planetesimals

def initialize_star_and_planetary_system(Mstar, Ndisk, Mdisk, Rmin, Rmax):

    converter=nbody_system.nbody_to_si(Mstar, Rmin)
    disk_massfraction = Mdisk/Mstar
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter,
                              densitypower=1.5, 
                              Rmin=Rmin.value_in(units.AU), 
                              Rmax=Rmax.value_in(units.AU),
                              q_out=1.0, discfraction=disk_massfraction).result
    #disk.h_smooth= Rmin/Ndisk

    star = Particles(1)
    star.mass = Mstar
    star.radius = 1| units.RSun
    star.position = (0,0,0) | units.AU
    star.velocity = (0,0,0) | units.kms
    planets = make_planets_oligarch.new_system(Mstar, star.radius,
                                               Rmin, Rmax, Mdisk)
    star.add_particles(planets[0].planets)
    print star
    
    return star, disk

def main(Mstar, Ndisk, fmdisk, Rmin, Rmax, t_end, n_steps):

    Mdisk = fmdisk * Mstar
    converter=nbody_system.nbody_to_si(Mstar, Rmax)
    star_and_planets, disk = initialize_star_and_planetary_system(Mstar, Ndisk, Mdisk, Rmin, Rmax)

    hydro=Fi(converter)
#    hydror = Fi(channel_type="ibis", hostname="galgewater")
 
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.timestep=0.125 | units.yr  

    hydro.gas_particles.add_particles(disk)

    gravity = Hermite(converter)
    gravity.particles.add_particles(star_and_planets)

    planet_attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "semimajor_axis", "eccentricity"]
    disk_attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "u", "rho", "h_smooth"]
    channel_to_planets = gravity.particles.new_channel_to(star_and_planets)
    channel_from_hydro_to_framework = hydro.particles.new_channel_to(disk,
                                                                     attributes=disk_attributes)
    channel_from_hydro_to_framework.copy()

    moving_bodies = ParticlesSuperset([star_and_planets, disk])
    moving_bodies.move_to_center()

    index = 0
    filename = "planetary_system_i{0:04}.amuse".format(index)
    write_set_to_file(star_and_planets, filename, 'amuse', attribute_names=planet_attributes,
                      append_to_file=False)
    write_set_to_file(disk, filename, 'amuse', attribute_names=disk_attributes)
    
    gravity_hydro = bridge.Bridge(use_threading=False)
    gravity_hydro.add_system(gravity, (hydro,) )
    gravity_hydro.add_system(hydro, (gravity,) )

    Etot_init = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)
    gravity_hydro.timestep = dt/10.
    while time < t_end:
        time += dt
        gravity_hydro.evolve_model(time)

        Etot_prev_se = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy

        channel_to_planets.copy()
        channel_from_hydro_to_framework.copy()

        for pi in star_and_planets[1:]:
            a, e = calculate_orbital_elements(star_and_planets[0], pi)
            pi.semimajor_axis = a
            pi.eccentricity = e
        
        index += 1
        filename = "planetary_system_i{0:04}.amuse".format(index)
        write_set_to_file(star_and_planets, filename, 'amuse', attribute_names=planet_attributes,
                          append_to_file=False)
        write_set_to_file(disk, filename, 'amuse', attribute_names=disk_attributes)

        Ekin = gravity_hydro.kinetic_energy 
        Epot = gravity_hydro.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

    gravity_hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 100,
                      help="number of diagnostics time steps [10]")
    result.add_option("--Ndisk", dest="Ndisk", type="int",default = 100000,
                      help="number of disk particles [100]")
    result.add_option("--Mstar", unit=units.MSun,
                      dest="Mstar", type="float",default = 1.73|units.MSun,
                      help="stellar mass [%default]")
    result.add_option("--fmdisk", dest="fmdisk", type="float",default = 0.01,
                      help="disk mass [%default]")
    result.add_option("--Rmin", unit=units.AU,
                      dest="Rmin", type="float",default = 1.0|units.AU,
                      help="minimal disk radius [%default]")
    result.add_option("--Rmax", unit=units.AU,
                      dest="Rmax", type="float",default = 100|units.AU,
                      help="maximal disk radius [%defualt]")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 1000.0|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


