import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from optparse import OptionParser
from amuse.couple import bridge

def main(Mstar=1, Ndisk=100, Mdisk= 0.001, Rmin=1, Rmax=100, t_end=10, n_steps=10, filename="nbody.hdf5"):
#    numpy.random.seed(111)
    Mstar = Mstar | units.MSun
    Mdisk = Mdisk | units.MSun
    Rmin = Rmin | units.AU
    Rmax = Rmax | units.AU
    t_end = t_end | units.yr

#    converter=nbody_system.nbody_to_si(Mdisk, Rmax)
    converter=nbody_system.nbody_to_si(Mstar, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, densitypower=1.5, 
                              Rmin=Rmin.value_in(units.AU), 
                              Rmax=Rmax.value_in(units.AU),q_out=1.0,
                              discfraction=1.0).result
 

    disk.move_to_center()

    center_of_mass = disk.center_of_mass()
    print(center_of_mass)
    a_bump = 50 
    Rbump = 10 | units.AU
    bump_location = center_of_mass + ([0, a_bump, 0] | units.AU)
    bump_particles = disk.select(lambda r: (center_of_mass-r).length()<Rbump,["position"])
    Mbump =bump_particles.mass.sum()

    Nbump = len(bump_particles)
    print("Nbump=", len(bump_particles), Mbump, Rbump, Nbump)
    print("Mass =", Mstar, Mdisk, disk.mass.sum().in_(units.MSun), disk.mass.sum()/Mstar)
    bump = new_plummer_gas_model(Nbump, convert_nbody=nbody_system.nbody_to_si(0.1*Mbump, Rbump))

    a_bump = a_bump | units.AU
    bump.x += a_bump
    r_bump = a_bump 
    inner_particles = disk.select(lambda r: (center_of_mass-r).length()<a_bump,["position"])
    M_inner = inner_particles.mass.sum() + Mstar

    v_circ = (constants.G*M_inner*(2./r_bump - 1./a_bump)).sqrt().value_in(units.kms)
    bump.velocity += [0, v_circ, 0] | units.kms
    print("MM=", M_inner, v_circ)

    disk.add_particles(bump)

    disk.h_smooth= Rmin/(Ndisk+Nbump)
    star=Particles(1)
    star.mass=Mstar
    star.radius=1. | units.RSun
    star.position = [0, 0, 0] | units.AU
    star.velocity = [0, 0, 0] | units.kms

#    star2=Particles(1)
#    star2.mass=Mstar
#    star2.radius=1. | units.RSun
#    star2.position = [1, 0, 0] | units.AU
#    star2.velocity = [0, 30, 0] | units.kms
 
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
    gravity.particles.add_particles(star)
#    gravity.particles.add_particles(star2)
#    gravity.particles.move_to_center()

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(star)
    channel_from_hydro_to_framework = hydro.particles.new_channel_to(disk)

    # make a repository of the star+disk particles in moving bodies
#    moving_bodies = star.union(disk)
    moving_bodies = ParticlesSuperset([star, disk])

    moving_bodies.move_to_center()

    write_set_to_file(moving_bodies, filename, 'hdf5')
    
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

        channel_from_gravity_to_framework.copy()
        channel_from_hydro_to_framework.copy()
        write_set_to_file(moving_bodies, filename, 'hdf5')

        Ekin = gravity_hydro.kinetic_energy 
        Epot = gravity_hydro.potential_energy
        Etot = Ekin + Epot
        print("T=", time, end=' ') 
        print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot) 
        print("star=", star)
        Etot_prev = Etot

    gravity_hydro.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 10,
                      help="number of diagnostics time steps [10]")
    result.add_option("-f", dest="filename", default = "gravhydro.hdf5",
                      help="output filename [gravhydro.hdf5]")
    result.add_option("--Ndisk", dest="Ndisk", type="int",default = 100,
                      help="number of disk particles [100]")
    result.add_option("--Mstar", dest="Mstar", type="float",default = 1,
                      help="stellar mass [1] MStar")
    result.add_option("--Mdisk", dest="Mdisk", type="float",default = 0.001,
                      help="disk mass [0.001] MStar")
    result.add_option("--Rmin", dest="Rmin", type="float",default = 1.0,
                      help="minimal disk radius [1] in AU")
    result.add_option("--Rmax", dest="Rmax", type="float",default = 100,
                      help="maximal disk radius [150] in AU")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] year")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

