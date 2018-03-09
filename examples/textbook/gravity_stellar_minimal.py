from amuse.lab import *
from amuse.units.optparse import OptionParser
    
def main(N, t_end, W0, Rvir, Mmin, Mmax, z):

###BOOKLISTSTART1###
    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter = nbody_system.nbody_to_si(Mtot_init, Rvir)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(convert_nbody=converter)
    
    stellar = SSE()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)
    bodies.radius = stellar.particles.radius

    gravity = ph4(converter)
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)

    channel_from_stellar_to_framework \
        = stellar.particles.new_channel_to(bodies)
    channel_from_stellar_to_gravity \
        = stellar.particles.new_channel_to(gravity.particles)
    channel_from_gravity_to_framework \
        = gravity.particles.new_channel_to(bodies)
###BOOKLISTSTOP1###

    filename = "bodies.hdf5"
    
###BOOKLISTSTART2###
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    dE_gr = 0 | Etot_init.unit
    time = 0.0 | t_end.unit
    dt = stellar.particles.time_step.amin()
    while time < t_end:

        dt = min(dt, t_end-time)
        
        stellar.evolve_model(time + dt/2)
        channel_from_stellar_to_gravity.copy()

        Etot_gr = gravity.kinetic_energy + gravity.potential_energy
        gravity.evolve_model(time + dt)
        dE_gr += (gravity.kinetic_energy + gravity.potential_energy - Etot_gr)

        stellar.evolve_model(time + dt)
        channel_from_stellar_to_gravity.copy()
        channel_from_gravity_to_framework.copy()

        time += dt
        write_set_to_file(bodies, filename)
###BOOKLISTSTOP2###

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot
    Mtot = bodies.mass.sum()
    print("T=", time, end=' ') 
    print("M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")", end=' ')
    print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
    print("dE/E=", (Etot_init-Etot)/Etot, end=' ')
    print("(dE[gr]/E=", dE_gr/Etot, ",", end=' ') 
    print("dE[se]/E=", (Etot_init-Etot-dE_gr)/Etot, ")")
    Etot_init -= dE

    gravity.stop()
    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [100]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mmax", type="float",default = 100 | units.MSun,
                      help="maximal stellar mass [100] MSun")
    result.add_option("-m", unit=units.MSun,
                      dest="Mmin", type="float",default = 0.1|units.MSun,
                      help="minimal stellar mass [0.1] MSun")
    result.add_option("-R", unit=units.parsec,
                      dest="Rvir", type="float",default = 1.0|units.parsec,
                      help="cluser virial radius [1] in parsec")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 1.0|units.Myr,
                      help="end time of the simulation [1] Myr")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="structure parameter for King model [7]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

