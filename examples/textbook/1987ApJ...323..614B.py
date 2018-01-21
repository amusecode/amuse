from amuse.lab import *
from amuse.io import store
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

def collide_two_stars(t_end, distance, offset, v_vesc, nsteps):
    filename = "Hydro_AM06MSun.h5"
    try:
        pstar = read_set_from_file(filename, format='hdf5')
    except:
        from local_star_to_sph import evolve_star_and_convert_to_sph
        mass = 0.6|units.MSun
        age = 8 | units.Gyr
        omega = 0|units.s**-1
        Nsph = 1000
        pstar, pcore = evolve_star_and_convert_to_sph(mass, age, omega, Nsph)

    print pstar
    pmass = pstar.mass.sum()
    print pmass.in_(units.MSun)

    filename = "Hydro_BM06MSun.h5"
    try:
        sstar = read_set_from_file(filename, format='hdf5')
    except:
        from local_star_to_sph import evolve_star_and_convert_to_sph
        mass = 0.6|units.MSun
        age = 8 | units.Gyr
        omega = 0|units.s**-1
        Nsph = 1000
        sstar, score = evolve_star_and_convert_to_sph(mass, age, omega, Nsph)
    print sstar
    smass = sstar.mass.sum()
    print smass.in_(units.MSun)

    import numpy
    v_esc = numpy.sqrt(2*constants.G*pmass/distance)
    velocity = v_vesc*v_esc
    
    sstar.x += distance
    sstar.y += offset
    sstar.vx -= velocity

    sph_particles = Particles(0)
    sph_particles.add_particles(pstar)
    sph_particles.add_particles(sstar)
    sph_particles.move_to_center()
    converter=nbody_system.nbody_to_si(pmass, distance)

    hydro = Fi(converter)

    dt = t_end/float(nsteps)

    hydro.gas_particles.add_particles(pstar)
    hydro.gas_particles.add_particles(sstar)
    Etot_init = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    to_framework = hydro.gas_particles.new_channel_to(sph_particles)

    output_file = "1987ApJ...323..614B.h5"
    write_set_to_file(sph_particles.savepoint(0.0 | t_end.unit), output_file, "hdf5", append_to_file=False)

    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt
        print time
        hydro.evolve_model(time)
        to_framework.copy()

        """
        from amuse.plot import sph_particles_plot
        sph_particles_plot(hydro.particles)
        from matplotlib import pyplot
        pyplot.show()
        """

        write_set_to_file(sph_particles.savepoint(time), output_file, "hdf5")

        """
        Ekin = hydro.kinetic_energy 
        Epot = hydro.potential_energy
        Eth = hydro.thermal_energy
        Etot = Ekin + Epot + Eth
        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), 
        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot
        """

    from amuse.plot import sph_particles_plot
    sph_particles_plot(hydro.particles)
    from matplotlib import pyplot
    pyplot.show()

    hydro.stop()

    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.hour,
                      dest="t_end", type="float", default = 10|units.hour,
                      help="end time of the simulation [%default]")
    result.add_option("-r", unit=units.MSun,
                      dest="distance", type="float", default = 1.4|units.RSun,
                      help="initial sepration [%default]")
    result.add_option("-o", unit=units.RSun,
                      dest="offset", type="float", default = 0|units.RSun,
                      help="offset sepration [%default]")
    result.add_option("-v", 
                      dest="v_vesc", type="float", default = 1.7,
                      help="impact velocity wrt escape speed [%default]")
    result.add_option("-n", 
                      dest="nsteps", type="int", default = 10,
                      help="number of steps [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    collide_two_stars(o.t_end, o.distance, o.offset, o.v_vesc, o.nsteps)

