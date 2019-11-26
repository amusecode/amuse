from __future__ import print_function
from amuse.lab import *
#from amuse.io import store
#from amuse.community.seba.interface import SeBa
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

def resolve_close_encounter(time, bodies, Johannes):
    Johannes.initialize_from_particles(bodies)
    rcom = bodies.center_of_mass()
    vcom = bodies.center_of_mass_velocity()

    a, e = Johannes.get_elements()
    p = Johannes.get_periastron()
    print("Close encounter at t=",  time.in_(units.Myr), "a=", a.in_(units.AU), "e=", e, "p=", p.in_(units.AU), "M=", bodies.mass.max().in_(units.MSun), bodies.mass.min().in_(units.MSun), "at d=", rcom.in_(units.parsec), "with v=", vcom.in_(units.kms))

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

    bodies[0].mass = bodies[0].stellar_mass + bodies[0].disk_mass + bodies[0].accreted_mass
    bodies[1].mass = bodies[1].stellar_mass + bodies[1].disk_mass + bodies[1].accreted_mass

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
    
def main(N, Rvir, Qvir, Fd):

    filename= 'Cl_N%g_R%gpc_Q%g_F%g.h5'%(N, Rvir.value_in(units.parsec), Qvir, Fd)
    t_end = 1.0 | units.Myr
    dt = 0.1 | units.Myr

    Mmax = 100 | units.MSun
    masses = new_kroupa_mass_distribution(N, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_fractal_cluster_model(N=N, fractal_dimension=Fd, 
                                           convert_nbody=converter)
    bodies.scale_to_standard(converter, virial_ratio=Qvir)
    bodies.stellar_mass = masses
    bodies.disk_mass = 0.1*bodies.stellar_mass
    bodies.mass = bodies.stellar_mass + bodies.disk_mass
    bodies.accreted_mass = 0 | units.MSun
    bodies.disk_radius = 400 | units.AU
    bodies.radius = 10 * bodies.disk_radius

    gravity = ph4(converter)
    gravity.parameters.epsilon_squared = (100|units.AU)**2
    gravity.particles.add_particles(bodies)

    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    Johannes = Kepler(converter)
    Johannes.initialize_code()

    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5', append_to_file=False)
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    Nenc = 0
    dEk_enc = zero    
    dEp_enc = zero
    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt

        gravity.evolve_model(time)
        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy
        while stopping_condition.is_set():
            channel_from_gd_to_framework.copy()
            Ek_enc = gravity.kinetic_energy 
            Ep_enc = gravity.potential_energy
            for ci in range(len(stopping_condition.particles(0))): 
                particles_in_encounter = Particles(particles=[stopping_condition.particles(0)[ci], stopping_condition.particles(1)[ci]])
                local_particles_in_encounter = particles_in_encounter.get_intersecting_subset_in(bodies)

                resolve_close_encounter(gravity.model_time, local_particles_in_encounter, Johannes)
                Nenc+=1
                print("At time=", gravity.model_time.value_in(units.Myr), "Nenc=", Nenc, "Rdisk=", local_particles_in_encounter.disk_radius.in_(units.AU))
                channel_from_framework_to_gd.copy_attributes(["radius"])
            dEk_enc += Ek_enc - gravity.kinetic_energy 
            dEp_enc += Ep_enc - gravity.potential_energy

            gravity.evolve_model(time)

        channel_from_framework_to_gd.copy_attributes(["mass"])

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_prev-Etot
        dE_se = Etot_prev_se-Etot
        Mtot = bodies.mass.sum()
        print("T=", time, end=' ') 
        print("M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")", end=' ')
        print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, end=' ') 
        print("(dE[SE]=", dE_se/Etot, ")")
        print("dE(enc)=", dEk_enc, dEp_enc)
        Etot_init -= dE
        Etot_prev = Etot

    gravity.stop()
    Johannes.stop()
    
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
    main(**o.__dict__)

