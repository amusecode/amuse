"""
   Simulate the hydrodynamial evolve a disk with a single bump around a star
"""
from amuse.lab import *
from amuse.io import store

def main(Mstar = 1|units.MSun,
         Ndisk=100, Mdisk=0.9|units.MSun, 
         Rmin=1.0|units.AU, Rmax=100.0|units.AU, 
         Mbump=0.1|units.MSun,Rbump=10.0|units.AU, abump=10|units.AU,
         t_end=1, n_steps=10):

    converter=nbody_system.nbody_to_si(Mdisk, Rmin)
    from amuse.ext.protodisk import ProtoPlanetaryDisk
    bodies = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, 
                                densitypower=1.5, Rmin=1.0, 
                                Rmax=Rmax/Rmin, q_out=1.0,
                                discfraction=1.0).result
    Mdisk = bodies.mass.sum()
    bodies.move_to_center()
    com = bodies.center_of_mass()

    mm = Mdisk/float(Ndisk)
    Nbump = Mbump/mm
    bump = new_plummer_gas_model(Nbump, convert_nbody=nbody_system.nbody_to_si(Mbump, Rbump))

    bump.x += abump
    r_bump = abump 
    inner_particles = bodies.select(lambda r: (com-r).length()<abump,["position"])
    M_inner = inner_particles.mass.sum() + Mstar
    v_circ = (constants.G*M_inner*(2./r_bump - 1./abump)).sqrt().value_in(units.kms)
    bump.velocity += [0, v_circ, 0] | units.kms
    bodies.add_particles(bump)

    star=Particles(1)
    star.mass=Mstar
    star.radius= Rmin
    star.position = [0, 0, 0] | units.AU
    star.velocity = [0, 0, 0] | units.kms

    import math
    P_bump = (abump**3*4*math.pi**2/(constants.G*(Mbump+Mstar))).sqrt()
    t_end *= P_bump

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(bodies)
    hydro.dm_particles.add_particles(star)
    Etot_init = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy

    particles = ParticlesSuperset([star, bodies])
    particles.move_to_center()
    particles.new_channel_to(hydro.particles).copy()
    bodies.h_smooth = Rmin # for the plotting routine
    channel_to_star = hydro.dm_particles.new_channel_to(star)
    channel_to_bodies = hydro.gas_particles.new_channel_to(bodies)

    write_set_to_file(star, "stars.hdf5","hdf5")
    write_set_to_file(bodies, "hydro.hdf5","hdf5")

    time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)
    while time < t_end:
        time += dt
        hydro.evolve_model(time)

        channel_to_star.copy()
        channel_to_bodies.copy()
        write_set_to_file(star, "stars.hdf5","hdf5")
        write_set_to_file(bodies, "hydro.hdf5","hdf5")
        star.radius = Rmin

        from hydro_sink_particles import hydro_sink_particles
        lost = hydro_sink_particles(star, bodies)
        if len(lost)>0:
            hydro.particles.remove_particles(lost)
            hydro.particles.synchronize_to(particles)
            print "Disk=", hydro.model_time, len(bodies), len(lost), lost.mass.sum(), star.mass

        Ekin = hydro.kinetic_energy 
        Epot = hydro.potential_energy
        Eth = hydro.thermal_energy
        Etot = Ekin + Epot + Eth
        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), 
        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot
        print "Star=", hydro.model_time, star[0].mass, star[0].position

    hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Ndisk", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-n", dest="n_steps", type="int",default = 10,
                      help="number of steps [10]")
    result.add_option("-t", 
                      dest="t_end", type="float", default = 1,
                      help="end time of the simulation in bump orbits")
    result.add_option("-M", dest="Mstar", type="float", default = 1|units.MSun,
                      help="Mass of the central star [%default]")
    result.add_option("--Mdisk", dest="Mdisk", type="float", 
                      default = 0.9|units.MSun,
                      help="Mass of the disk [%default]")
    result.add_option("-r", unit=units.AU,
                      dest="Rmin", type="float", default = 10 |units.AU,
                      help="inner disk radius [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rmax", type="float", default = 100 | units.AU,
                      help="outer disk radius [%default]")
    result.add_option("--Mbump", unit=units.MSun,
                      dest="Mbump", type="float", default = 0.1 | units.MSun,
                      help="bump mass [%default]")
    result.add_option("--Rbump", unit=units.AU,
                      dest="Rbump", type="float", default = 10 | units.AU,
                      help="bump radius [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="abump", type="float", default = 20 | units.AU,
                      help="distance of bump from star [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

