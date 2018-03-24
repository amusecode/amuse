import numpy
import cPickle
from amuse.lab import *
from amuse.community.fastkick.interface import FastKick
from amuse.ext.relax_sph import relax
from amuse.ext.spherical_model import new_gas_plummer_distribution
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

###BOOKLISTSTART1###
def check_energy_conservation(system, i_step, time, n_steps):
    unit = units.J
    U = system.potential_energy.value_in(unit)
    Q = system.thermal_energy.value_in(unit)
    K = system.kinetic_energy.value_in(unit)
    print "Step {0} of {1}, t={2}: U={3:.2e}, Q={4:.2e}, K={5:.2e} {5}".format(
        i_step, n_steps, time.as_quantity_in(units.yr), U, Q, K, unit)
###BOOKLISTSTOP1###

from amuse.couple.bridge import Bridge
def local_relax(gas_particles, hydro, gravity_field=None,
                monitor_func=check_energy_conservation,
                bridge_options=dict()):
    """
    Relax a set of SPH particles by evolving it with a hydrodynamics code, while 
    imposing critical damping on the particle velocities.
    
    :argument gas_particles:  The set of SPH particles
    :argument hydro:          The hydrodynamics code
    :argument gravity_field   Background gravitational field, must support get_gravity_at_point
    :argument monitor_func    For monitoring progress each step. User-defined function or "energy"
    :argument bridge_options: Keyword options passed to Bridge
    """
    if monitor_func == "energy":
        monitor_func = monitor_energy
    t_end_in_t_dyn = 0.1 # Relax for this many dynamical timescales
    t_end = t_end_in_t_dyn \
              * gas_particles.dynamical_timescale(mass_fraction=0.9)
    n_steps = 10
    velocity_damp_factor = 1.0 - (2.0*numpy.pi*t_end_in_t_dyn) \
                                   /n_steps # Critical damping
    
    in_hydro = hydro.gas_particles.add_particles(gas_particles)
    if gravity_field is None:
        system = hydro
    else:
        system = Bridge(timestep=(t_end/n_steps).as_quantity_in(units.yr),
                        **bridge_options)
        system.add_system(hydro, [gravity_field])
    
    for i_step, time in enumerate(t_end * numpy.linspace(1.0/n_steps,
                                                         1.0, n_steps)):
        system.evolve_model(time)
        hydro.gas_particles.velocity = velocity_damp_factor \
                                         * hydro.gas_particles.velocity
        monitor_func(system, i_step, time, n_steps)
    
    return in_hydro.copy()

###BOOKLISTSTART2###
def relax_gas_and_stars(stars, gas):
    dynamical_timescale = gas.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.parsec)
    hydro = Fi(converter, mode='openmp', redirection="file",
               redirect_file="fi.log")
    hydro.parameters.timestep = dynamical_timescale / 100
    hydro.parameters.eps_is_h_flag = True
    
    gravity_field_code = FastKick(converter, mode="cpu", number_of_workers=2)
    gravity_field_code.parameters.epsilon_squared = (0.01 | units.parsec)**2
    gravity_field_code.particles.add_particles(stars)
    monitoring_function = check_energy_conservation
    relaxed_gas = relax(gas, hydro, gravity_field=gravity_field_code, 
                        monitor_func=monitoring_function,
                        bridge_options=dict(verbose=True, use_threading=False))
    gravity_field_code.stop()
    return hydro
###BOOKLISTSTOP2###

def generate_initial_conditions(
        number_of_stars = 100,
        number_of_gas_particles = 10**5,
        star_formation_efficiency = 0.1,
        virial_radius = 0.33 | units.parsec,
        virial_ratio = 1.0,
        use_fractal = False):
    
    numpy.random.seed(12345678)
    seed_fractal = 312357271
    
    masses = new_salpeter_mass_distribution(number_of_stars,
                                            mass_min=1|units.MSun)
#    masses = new_kroupa_mass_distribution(number_of_stars
    total_stellar_mass = masses.sum()
    total_mass = total_stellar_mass / star_formation_efficiency
    converter = nbody_system.nbody_to_si(total_mass, virial_radius)
    if use_fractal:
        stars = new_fractal_cluster_model(number_of_stars,
                                          convert_nbody=converter,
                                          do_scale=False,
                                          fractal_dimension=1.6,
                                          random_seed=seed_fractal)
    else:
        stars = new_plummer_model(number_of_stars, convert_nbody=converter,
                                  do_scale=False)
    stars.mass = masses
    stars.move_to_center()
    print "scaling positions to match virial_radius"
    stars.position *= virial_radius / stars.virial_radius()
    print "scaling velocities to match virial_ratio"
    stars.velocity *= numpy.sqrt(virial_ratio * converter.to_si(0.5|nbody_system.energy) * star_formation_efficiency / stars.kinetic_energy())
    
    print "new_gas_plummer_distribution"
    gas = new_gas_plummer_distribution(
        number_of_gas_particles, 
        total_mass = (total_mass - total_stellar_mass), 
        virial_radius = virial_radius, 
        type = "fcc")
    gas.h_smooth = 0.0 | units.parsec
    # eat away gas.
    mgas = gas[0].mass
    print "Ngas=", len(gas)
    for si in stars:
        m = si.mass
        nremoved = 0
        while m>mgas:
            gi = si.as_set().nearest_neighbour(gas)
            nremoved += 1
            m-= gi.mass
            gas -= gi
        print "removed:", nremoved, si.mass.in_(units.MSun)
    print "Ngas=", len(gas)
        
    filename = "YSC_{0}_stars{1}_gas{2}k_" \
                 .format("fractal" if use_fractal else "plummer",
        number_of_stars, number_of_gas_particles/1000)
    print "Writing initial conditions to", filename, "+ stars/gas.amuse"
    write_set_to_file(stars, filename+"stars.amuse", "amuse",
                      append_to_file=False)
    write_set_to_file(gas, filename+"gas.amuse", "amuse", append_to_file=False)
    with open(filename+"info.pkl", "wb") as outfile:
        cPickle.dump([converter], outfile)
    return stars, gas, filename             

def make_map(sph,N=100,L=1):
    x,y=numpy.indices( ( N+1,N+1 ))
    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.parsec(x)
    y=units.parsec(y)
    z=units.parsec(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))
    return rho

def plot_hydro_and_stars(hydro, stars):
    x_label = "x [pc]"
    y_label = "y [pc]"
    from prepare_figure import *
    fig = single_frame(x_label, y_label, logx=False, logy=False,
                       xsize=12, ysize=12)

    L = 6
    rho=make_map(hydro,N=200,L=L)
    cax = pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
                        extent=[-L/2,L/2,-L/2,L/2],vmin=1,vmax=5,
                        origin="lower")
    cbar = fig.colorbar(cax, orientation='vertical', fraction=0.045)
    cbar.set_ticks([1.0, 1.5, 2.0, 4.0, 4.5, 5.0])
    cbar.set_label('$\log_{10} (\mathrm{n}) \,\,\, [\mathrm{amu/cm}^3]$',
                   rotation=270)
    
    cm = pyplot.cm.get_cmap('RdBu')
    m = 50*numpy.log10(stars.mass/stars.mass.min())
    c = numpy.sqrt(stars.mass/stars.mass.max())
    pyplot.scatter(-stars.x.value_in(units.parsec),
                   -stars.y.value_in(units.parsec),
                   c=c, s=m, lw=0, cmap=cm)

    save_file = 'plot_relaxed_gas_and_star.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
    
if __name__ == "__main__":
    stars, gas, filename = generate_initial_conditions()
    hydro = relax_gas_and_stars(stars, gas)
    plot_hydro_and_stars(hydro, stars)
    hydro.stop()


    
