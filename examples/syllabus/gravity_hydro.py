"""
   example code for bridging a gravity solver with a hydrodynamics solver
"""
import math
import numpy
from amuse.lab import *
from amuse.couple import bridge

def print_diagnostics(model_time, bh, star, converter):
    kep = Kepler(converter)
    kep.initialize_code()

    mass = bh.mass + star.mass
    pos = bh.position - star.position
    vel = bh.velocity - star.velocity
    print mass, pos, vel
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2], 
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    r = kep.get_separation()
    E,J = kep.get_integrals()
    kep.stop()
    print "Orbital elements for star ", model_time, a, e, bh.mass, star.mass
    return a, e

def orbital_period(a, Mtot) :
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()

def main(a=1|units.AU, ecc=0, t_end=1.|units.yr, n_steps=100, 
         Rgas=1|units.AU, Mgas=1|units.MSun, Ngas=1000):
    stars = Particles(2)
    stars.mass = [3.2, 3.1] | units.MSun
    Porb = orbital_period(a, stars.mass.sum())
    stars[0].position = (0,0,0) | units.AU
    stars[0].velocity = (0,0,0) | units.kms
    apo = a*(1+ecc)
    vc = math.sqrt((1-ecc)/(1+ecc)) * (constants.G*stars.mass.sum()/apo).sqrt()
    stars[1].position = (a.value_in(units.AU),0,0) | units.AU
    stars[1].velocity = (0,vc.value_in(units.kms),0) | units.kms
    stars.move_to_center()
    converter=nbody_system.nbody_to_si(stars.mass.sum(), a)
    gravity = ph4(converter, redirection="none")
    gravity.particles.add_particles(stars)
    gravity.parameters.epsilon_squared = (10|units.RSun)**2
    Ed0_tot = gravity.kinetic_energy + gravity.potential_energy

    channel_from_gravity = gravity.particles.new_channel_to(stars)
    channel_from_to_gravity = stars.new_channel_to(gravity.particles)

    dt = t_end/float(n_steps)

    converter=nbody_system.nbody_to_si(1.0|units.MSun, Rgas)
    ism = new_plummer_gas_model(Ngas, convert_nbody=converter)
    ism.move_to_center()
    ism = ism.select(lambda r: r.length()<2*a,["position"])
    print "Ngas=", len(ism)

    hydro = Fi(converter, redirection="none")
    hydro.parameters.timestep = dt/8.
    hydro.parameters.epsilon_squared = (10|units.RSun)**2
#    hydro.parameters.epsilon_squared = 0.25*(float(Ngas))**(-0.666667) * a**2
    hydro.gas_particles.add_particles(ism)
    Eh0_tot = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    hydro.parameters.periodic_box_size = 10*Rgas

    channel_from_hydro = hydro.gas_particles.new_channel_to(ism)
    channel_from_to_hydro = ism.new_channel_to(hydro.gas_particles)

    model_time = 0 | units.Myr
    filename = "gravhydro.hdf5"
    write_set_to_file(stars.savepoint(model_time), filename, 'hdf5')
    write_set_to_file(ism, filename, 'hdf5')
#    write_set_to_file(ism.savepoint(model_time), filename, 'hdf5')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,) )
    gravhydro.add_system(hydro, (gravity,) )
    print "dt=", 0.01*dt, Porb/1024.
    gravhydro.timestep = min(dt, 2*hydro.parameters.timestep)

    snapshot_id = 0
    while model_time < t_end:
        snapshot_id += 1
        model_time += dt

        a, e = print_diagnostics(gravity.model_time, stars[0], stars[1], converter) 
        print "time=", model_time, a, e
        gravhydro.evolve_model(model_time)

        Ed_tot = gravity.kinetic_energy + gravity.potential_energy
        Eh_tot = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
        print "Energies:", Ed_tot/Ed0_tot, Eh_tot/Eh0_tot

        channel_from_gravity.copy()
        channel_from_hydro.copy()

        if snapshot_id%10==0:
            write_set_to_file(stars.savepoint(model_time), filename, 'hdf5')
            write_set_to_file(ism, filename, 'hdf5')
#        write_set_to_file(ism.savepoint(model_time), filename, 'hdf5')

    gravity.stop()
    hydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1000,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ngas", type="int", default = 1024,
                      help="number of gas particles [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mgas", type="float", default = 1|units.MSun,
                      help="Mass of the gas [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rgas", type="float", default = 1|units.AU,
                      help="Size of the gas distribution [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="a", type="float", default = 0.1|units.AU,
                      help="initial orbital separation [%default]")
    result.add_option("-e", dest="ecc", type="float", default = 0.0,
                      help="initial orbital eccentricity [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 10|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
