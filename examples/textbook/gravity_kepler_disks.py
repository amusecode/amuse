from amuse.lab import *
#from amuse.io import store
#from amuse.community.seba.interface import SeBa
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

###BOOKLISTSTART3###
def resolve_close_encounter(time, bodies):
    orbital_elements = orbital_elements_from_binary(bodies, G=constants.G)
    a = orbital_elements[2]
    e = orbital_elements[3]
    p = a*(1-e)
    print("Close encounter at t=",  time.in_(units.Myr), \
          "a=", a.in_(units.AU), "e=", e, "p=", p.in_(units.AU), \
          "M=", bodies.mass.max().in_(units.MSun), \
          bodies.mass.min().in_(units.MSun)) 
    truncate_disks_due_to_encounter(bodies, p)

def truncate_disks_due_to_encounter(bodies, p):
    q = bodies[1].mass/bodies[0].mass
    rtr_prim = 0.28*p / q**(0.32)
    rtr_sec = 0.28*p * q**(0.32)
    dm0 = truncate_disk_due_to_encounter(bodies[0], rtr_prim)
    dm1 = truncate_disk_due_to_encounter(bodies[1], rtr_sec)

    mtot = bodies.mass.sum()
    bodies[0].accreted_mass += dm1 * bodies[0].mass/mtot
    bodies[1].accreted_mass += dm0 * bodies[1].mass/mtot

    bodies[0].radius = min(bodies[0].radius, 0.5*p)
    bodies[1].radius = min(bodies[1].radius, 0.5*p)

    bodies[0].mass = bodies[0].stellar_mass + bodies[0].disk_mass \
                       + bodies[0].accreted_mass
    bodies[1].mass = bodies[1].stellar_mass + bodies[1].disk_mass \
                       + bodies[1].accreted_mass
###BOOKLISTSTOP3###

def stripped_disk_mass(body, dr):
    rold = body.disk_radius
    rnew = rold - dr
    dm = body.disk_mass * (rold**0.5-rnew**0.5)/rold**0.5
    return max(0|units.MSun, dm)

def truncate_disk_due_to_encounter(body, r_tr):
    dr = max(0|units.AU, body.disk_radius-r_tr)
    dm = stripped_disk_mass(body, dr)
    body.disk_radius -= dr
    body.disk_mass -= dm
    return dm

###BOOKLISTSTART2###
def evolve_system_to(time, gravity, bodies, stopping_condition,
                     channel_from_gravity, channel_to_gravity,
                     energy_tolerance = 1.e-10 | units.erg):
                     
    gravity.evolve_model(time)

    while stopping_condition.is_set():
        channel_from_gravity.copy()
        Ek_enc = gravity.kinetic_energy 
        Ep_enc = gravity.potential_energy
        for ci in range(len(stopping_condition.particles(0))): 
            bodies_in_enc \
                = Particles(particles=[stopping_condition.particles(0)[ci],
                                       stopping_condition.particles(1)[ci]])
            local_bodies_in_enc \
                = bodies_in_enc.get_intersecting_subset_in(bodies)

            resolve_close_encounter(gravity.model_time, local_bodies_in_enc)
            print("At time=", gravity.model_time.value_in(units.Myr), \
                "Rdisk=", local_bodies_in_enc.disk_radius.in_(units.AU))
            channel_to_gravity.copy_attributes(["radius"])
            assert abs(Ek_enc - gravity.kinetic_energy) < energy_tolerance
            assert abs(Ep_enc - gravity.potential_energy) < energy_tolerance
        gravity.evolve_model(time)
    channel_to_gravity.copy_attributes(["mass"])
###BOOKLISTSTOP2###

###BOOKLISTSTART1###
def main(N, Rvir, Qvir, Fd, t_end, filename):
    masses = new_kroupa_mass_distribution(N, 100|units.MSun)
    converter=nbody_system.nbody_to_si(masses.sum(),Rvir)
    bodies = new_fractal_cluster_model(N=N, fractal_dimension=Fd, 
                                       convert_nbody=converter)
    bodies.scale_to_standard(converter, virial_ratio=Qvir)
    bodies.stellar_mass = masses
    bodies.disk_mass = 0.01*bodies.stellar_mass
    bodies.mass = bodies.stellar_mass + bodies.disk_mass
    bodies.accreted_mass = 0 | units.MSun
    bodies.disk_radius = 400 | units.AU
    bodies.radius = 10 * bodies.disk_radius

    gravity = ph4(converter, number_of_workers=2)
    gravity.parameters.epsilon_squared = (100|units.AU)**2
    gravity.particles.add_particles(bodies)
    channel_from_gravity = gravity.particles.new_channel_to(bodies)
    channel_to_gravity = bodies.new_channel_to(gravity.particles)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5',
                      append_to_file=False)
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init
###BOOKLISTSTOP1###

###BOOKLISTSTART0###
    dt = t_end/10.
    time = 0 | units.yr
    while gravity.model_time < t_end:
        time += dt
        evolve_system_to(time, gravity, bodies, stopping_condition,
                         channel_from_gravity, channel_to_gravity)
        write_set_to_file(bodies.savepoint(gravity.model_time),
                          filename, 'hdf5')
        Etot = gravity.kinetic_energy + gravity.potential_energy
        print("T=", gravity.model_time, end=' ') 
        print("E= ", Etot, "Q= ", \
              gravity.kinetic_energy/gravity.potential_energy)
        print("dE=", (Etot-Etot_init)/Etot, "ddE=", (Etot-Etot_prev)/Etot)
        Etot_init -= (Etot_prev-Etot)
        Etot_prev = Etot

    gravity.stop()
###BOOKLISTSTOP0###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 2000,
                      help="number of stars [%default]")
    result.add_option("-R", dest="Rvir", type="float",
                      unit=units.parsec, default = 0.5|units.parsec,
                      help="cluser virial radius [%default]")
    result.add_option("-Q", dest="Qvir", type="float",default = 0.5,
                      help="virial ratio [%default]")
    result.add_option("-F", dest="Fd", type="float",default = 1.6,
                      help="fractal dimension [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    t_end = 1|units.Myr
    filename= 'Cl_N%g_R%gpc_Q%g_F%g.h5' \
               % (o.N, o.Rvir.value_in(units.parsec), o.Qvir, o.Fd)
    main(o.N, o.Rvir, o.Qvir, o.Fd, t_end, filename)

