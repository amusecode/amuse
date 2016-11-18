"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
from amuse.lab import *
from optparse import OptionParser
    
def main(N=10, W0=7.0, t_end=10, Rvir=1, Mmin=0.1, Mmax= 100, z=0.02):
    t_end = t_end | units.Myr
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses

    stellar = SSE()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)
    stellar.commit_particles()

    bodies.radius = stellar.particles.radius

    gravity = ph4(converter)
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_se_to_framework.copy_attributes(["mass"])
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)

    bodies.scale_to_standard(convert_nbody=converter)
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    dE_gr = 0 | Etot_init.unit
    dE_se = 0 | Etot_init.unit
    time = 0.0 | t_end.unit
    dt = stellar.particles.time_step.amin()
    while time < t_end:
        time = min(time+dt, t_end)

        Etot_gr = gravity.kinetic_energy + gravity.potential_energy
        gravity.evolve_model(time)
        dE_gr += (gravity.kinetic_energy + gravity.potential_energy-Etot_gr)

        stellar.evolve_model(time)
        Etot_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["mass"])
        channel_from_framework_to_gd.copy_attributes(["mass"])

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot
    Mtot = bodies.mass.sum()
    print "T=", time, 
    print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
    print "E= ", Etot, "Q= ", Ekin/Epot,
    print "dE/E=", (Etot_init-Etot)/Etot,
    print "(dE[gr]/E=", dE_gr/Etot, ",", 
    print "dE[se]/E=", (Etot_init-Etot-dE_gr)/Etot, ")"
    Etot_init -= dE

    gravity.stop()
    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [100]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [100] MSun")
    result.add_option("-m", dest="Mmin", type="float",default = 0.1,
                      help="minimal stellar mass [0.1] MSun")
    result.add_option("-R", dest="Rvir", type="float",default = 1.0,
                      help="cluser virial radius [1] in parsec")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] Myr")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

