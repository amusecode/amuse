"""
   example code bridging a gravity solver with a hydrodynamics solver
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge

def evolve_binary_in_common_envelope(stars, envelope, t_end):
    R = stars.position.length()
    converter=nbody_system.nbody_to_si(stars.mass.sum(), R)

    gravity = ph4(converter)
    gravity.particles.add_particles(stars)

    channel_from_gravity = gravity.particles.new_channel_to(stars)
    channel_from_to_gravity = stars.new_channel_to(gravity.particles)

    n_steps = 10
    dt = t_end/float(n_steps)

    hydro = Fi(converter, redirection="none")
    tdyn = numpy.sqrt((0.05*R)**3/(constants.G*stars.mass.sum()))
    print "tdyn=", tdyn
    hydro.parameters.timestep = tdyn
    hydro.parameters.epsilon_squared = (1|units.RSun)**2
    hydro.gas_particles.add_particles(envelope)
    hydro.parameters.periodic_box_size = 100*R

    channel_from_hydro = hydro.gas_particles.new_channel_to(envelope)
    channel_from_to_hydro = envelope.new_channel_to(hydro.gas_particles)

    model_time = 0 | units.Myr
    filename = "XiTau_Hydro.amuse"
    write_set_to_file(stars.savepoint(model_time), filename, 'amuse',
                      append_to_file=False)
    write_set_to_file(envelope, filename, 'amuse')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,) )
    gravhydro.add_system(hydro, (gravity,) )
    gravhydro.timestep = min(dt, 10*hydro.parameters.timestep)

    while model_time < t_end:
        model_time += dt
        print "Time=", model_time.in_(units.day)
        gravhydro.evolve_model(model_time)

        channel_from_gravity.copy()
        channel_from_hydro.copy()

        write_set_to_file(stars.savepoint(model_time), filename, 'amuse')
        write_set_to_file(envelope, filename, 'amuse')
    gravity.stop()
    hydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.Myr, 
                      dest="t_end", type="float", default = 20|units.day,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    a = 0.133 | units.AU
    e = 0.0
    m1 = 3.2|units.MSun
    m2 = 3.1|units.MSun
    from amuse.ext.orbital_elements import new_binary_from_orbital_elements
    inner_binary = new_binary_from_orbital_elements(m1, m2, a, e,
                                                    G=constants.G)

    XiTau = read_set_from_file("Hydro_PrimaryStar_XiTau.amuse", "amuse")
    for ci in XiTau.history:
        if len(ci) == 1:
            XiTau_core = ci.copy()
        else:
            XiTau_envelope = ci.copy()

    ao = 1.23 | units.AU
    eo = 0.15
    m12 = inner_binary.mass.sum()
    m3 = XiTau_core.mass + XiTau_envelope.mass.sum()
    print "M3=", m3.in_(units.MSun)
    outer_binary = new_binary_from_orbital_elements(m12, m3, ao, eo,
                                                    G=constants.G)

    inner_binary.position += outer_binary[0].position
    inner_binary.velocity += outer_binary[0].velocity

    XiTau_envelope.position += outer_binary[1].position
    XiTau_envelope.velocity += outer_binary[1].velocity
    outer_binary[1].mass = XiTau_core.mass

    triple = Particles()
    triple.add_particles(inner_binary)
    triple.add_particle(outer_binary[1])

    print triple
    
    evolve_binary_in_common_envelope(triple, XiTau_envelope, o.t_end)

