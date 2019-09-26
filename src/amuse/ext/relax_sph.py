import numpy
from amuse.units import units, nbody_system
from amuse.couple.bridge import Bridge

from amuse.plot import native_plot, semilogy, loglog, xlabel, ylabel

def no_monitoring(system, i_step, time, n_steps):
    pass

def monitor_energy(system, i_step, time, n_steps):
    unit = units.J
    U = system.potential_energy.value_in(unit)
    Q = system.thermal_energy.value_in(unit)
    K = system.kinetic_energy.value_in(unit)
    print("Step {0}, t={1}: U={2:.2e}, Q={3:.2e}, K={4:.2e} {5}".format(
        i_step, time.as_quantity_in(units.yr), U, Q, K, unit))

class Memory:
    pass

def monitor_density_profile(system, i_step, time, n_steps, memory=Memory()):
    if i_step == 0:
        memory.xlimits = (None, None)
        memory.ylimits = (None, None)
    position = system.gas_particles.position - system.gas_particles.center_of_mass()
    loglog(position.lengths_squared(), system.gas_particles.density, 'gs')
    native_plot.title("{0}: t={1}".format(i_step, time.as_quantity_in(units.yr)))
    native_plot.xlim(memory.xlimits)
    native_plot.ylim(memory.ylimits)
    native_plot.pause(0.0001)
    memory.xlimits = native_plot.gca().get_xlim()
    memory.ylimits = native_plot.gca().get_ylim()
    if i_step == n_steps-1:
        native_plot.show(block=True)
    native_plot.cla()

def relax(gas_particles, hydro, gravity_field=None, monitor_func=no_monitoring,
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
    t_end_in_t_dyn = 2.5 # Relax for this many dynamical timescales
    t_end = t_end_in_t_dyn * gas_particles.dynamical_timescale(mass_fraction=0.9)
    n_steps = 250
    velocity_damp_factor = 1.0 - (2.0*numpy.pi*t_end_in_t_dyn)/n_steps # Critical damping
    
    in_hydro = hydro.gas_particles.add_particles(gas_particles)
    if gravity_field is None:
        system = hydro
    else:
        system = Bridge(timestep=(t_end/n_steps).as_quantity_in(units.yr), **bridge_options)
        system.add_system(hydro, [gravity_field])
    
    for i_step, time in enumerate(t_end * numpy.linspace(1.0/n_steps, 1.0, n_steps)):
        system.evolve_model(time)
        hydro.gas_particles.velocity = velocity_damp_factor * hydro.gas_particles.velocity
        monitor_func(system, i_step, time, n_steps)
    
    return in_hydro.copy()


if __name__ == "__main__":
    from amuse.io import write_set_to_file
    from amuse.ext.spherical_model import new_gas_plummer_distribution, new_plummer_distribution
    from amuse.community.gadget2.interface import Gadget2
    from amuse.community.fastkick.interface import FastKick
    gas = new_gas_plummer_distribution(1000, virial_radius=1|units.parsec, total_mass=1000|units.MSun, type="fcc")
    stars = new_plummer_distribution(10, virial_radius=1|units.parsec, total_mass=100|units.MSun, type="sobol")
    
    dynamical_timescale = gas.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.parsec)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.parameters.time_max = 3 * dynamical_timescale
    hydro.parameters.max_size_timestep = dynamical_timescale / 100
    hydro.parameters.time_limit_cpu = 1.0 | units.Gyr
    
    gravity_field_code = FastKick(converter)
    gravity_field_code.particles.add_particles(stars)
    relaxed_gas = relax(gas, hydro, gravity_field=gravity_field_code, monitor_func="energy", bridge_options=dict(verbose=True))
    gravity_field_code.stop()
    hydro.stop()
    write_set_to_file(relaxed_gas, "gas_relaxed.amuse", "amuse")
