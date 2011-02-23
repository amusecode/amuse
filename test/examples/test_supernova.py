import os.path
from amuse.test.amusetest import get_path_to_results
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles, Grid
from amuse.support.units import units, generic_unit_system, nbody_system, constants
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.exceptions import AmuseException
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.ext.star_to_sph import *

def inject_supernova_energy(gas_particles):
    inner = gas_particles.select(lambda pos : pos.length_squared() < 100.0 | units.RSun**2, ["position"])
    print len(inner), "innermost particles selected."
    print "Adding", (1.0e51 | units.erg) / inner.total_mass(), "of supernova " \
        "(specific internal) energy to each of them."
    inner.u += (1.0e51 | units.erg) / inner.total_mass()

def setup_stellar_evolution_model():
    out_pickle_file = os.path.join(get_path_to_results(), "super_giant_stellar_structure.pkl")
    if os.path.exists(out_pickle_file):
        return out_pickle_file
    
    stellar_evolution = MESA(redirection = "none")
    stars =  Particles(1)
    stars.mass = 10.0 | units.MSun
    stellar_evolution.initialize_module_with_default_parameters() 
    stellar_evolution.particles.add_particles(stars)
    stellar_evolution.initialize_stars()

    print "Evolving a MESA star with mass:", stellar_evolution.particles[0].mass
    try:
        while True:
            stellar_evolution.evolve_model()
    except AmuseException as ex:
        print "Evolved star to", stellar_evolution.particles[0].age
        print "Radius:", stellar_evolution.particles[0].radius
    
    pickle_stellar_model(stellar_evolution.particles[0], out_pickle_file)
    stellar_evolution.stop()
    return out_pickle_file

def run_supernova():
    # options:
    use_hydro_code = Gadget2 # Fi -or- Gadget2
    hydro_code_options = dict(number_of_workers=3) # e.g. dict(use_gl = True) or dict(redirection = "none")
    number_of_sph_particles = 30000
    t_end = 1.0e5 | units.s
    
    pickle_file = setup_stellar_evolution_model()
    
    print "Creating initial conditions from a MESA stellar evolution model..."
    core, gas_without_core, core_radius = convert_stellar_model_to_SPH(
        None, 
        number_of_sph_particles, 
        seed = 12345,
        pickle_file = pickle_file,
        with_core_particle = True
    )
    if len(core):
        print "Created", len(gas_without_core), "SPH particles and one 'core-particle':\n", core
        print "Setting gravitational smoothing to:", core_radius
    else:
        print "Warning: Only SPH particles created."
    
    inject_supernova_energy(gas_without_core)
    
    print "\nEvolving (SPH) to:", t_end
    n_steps = 100
    
    unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, constants.G, t_end)
    hydro_code = use_hydro_code(unit_converter, **hydro_code_options)
    
    try:
        hydro_code.parameters.timestep = t_end / n_steps
    except Exception as exc:
        if not "parameter is read-only" in str(exc): raise
    
    hydro_code.parameters.epsilon_squared = core_radius**2
    hydro_code.parameters.n_smooth_tol = 0.01 | units.none
    hydro_code.gas_particles.add_particles(gas_without_core)
    hydro_code.dm_particles.add_particles(core)
    
    times = [0.0] | units.s
    potential_energies = hydro_code.potential_energy.as_quantity_in(units.J).as_vector_with_length(1)
    kinetic_energies =   hydro_code.kinetic_energy.as_quantity_in(units.J).as_vector_with_length(1)
    thermal_energies =   hydro_code.thermal_energy.as_quantity_in(units.J).as_vector_with_length(1)
    for time, i_step in [(i*t_end/n_steps, i) for i in range(1, n_steps+1)]:
        hydro_code.evolve_model(time)
        times.append(time)
        potential_energies.append( hydro_code.potential_energy)
        kinetic_energies.append(   hydro_code.kinetic_energy)
        thermal_energies.append(   hydro_code.thermal_energy)
        hydro_plot(
            [-1.0, 1.0, -1.0, 1.0] * (350 | units.RSun),
            hydro_code,
            (100, 100),
            os.path.join(get_path_to_results(), "supernova_hydro_image{0:=03}.png".format(i_step))
        )
    
    energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
        os.path.join(get_path_to_results(), "supernova_energy_evolution.png"))

    hydro_code.stop()
    print "All done!\n"
    

def energy_plot(time, E_kin, E_pot, E_therm, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    plot(time, E_kin.as_quantity_in(units.erg), label='E_kin')
    plot(time, E_pot, label='E_pot')
    plot(time, E_therm, label='E_therm')
    plot(time, E_kin+E_pot+E_therm, label='E_total')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(loc=3)
    pyplot.savefig(figname)
    print "\nPlot of energy evolution was saved to: ", figname
    pyplot.close()

def hydro_plot(view, hydro_code, image_size, figname):
    """
    view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    hydro_code: hydrodynamics code in which the gas to be plotted is defined
    image_size: size of the output image in pixels (x, y)
    """
    if not HAS_MATPLOTLIB:
        return
    shape = (image_size[0], image_size[1], 1)
    size = image_size[0] * image_size[1]
    axis_lengths = [0.0, 0.0, 0.0] | units.m
    axis_lengths[0] = view[1] - view[0]
    axis_lengths[1] = view[3] - view[2]
    grid = Grid.create(shape, axis_lengths)
    grid.x += view[0]
    grid.y += view[2]
    speed = grid.z.reshape(size) * (0 | 1/units.s)
    rho, rhovx, rhovy, rhovz, rhoe = hydro_code.get_hydro_state_at_point(grid.x.reshape(size), 
        grid.y.reshape(size), grid.z.reshape(size), speed, speed, speed)
    
    min_v =  800.0 | units.km / units.s
    max_v = 3000.0 | units.km / units.s
    min_rho = 3.0e-9 | units.g / units.cm**3
    max_rho = 1.0e-5 | units.g / units.cm**3
    min_E = 1.0e11 | units.J / units.kg
    max_E = 1.0e13 | units.J / units.kg
    
    v_sqr = (rhovx**2 + rhovy**2 + rhovz**2) / rho**2
    E = rhoe / rho
    log_v = numpy.log((v_sqr / min_v**2).value_in(units.none)) / numpy.log((max_v**2 / min_v**2).value_in(units.none))
    log_rho = numpy.log((rho / min_rho).value_in(units.none)) / numpy.log((max_rho / min_rho).value_in(units.none))
    log_E = numpy.log((E / min_E).value_in(units.none)) / numpy.log((max_E / min_E).value_in(units.none))
    
    red   = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_rho)).reshape(shape)
    green = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_v)).reshape(shape)
    blue  = numpy.minimum(numpy.ones_like(rho.number), numpy.maximum(numpy.zeros_like(rho.number), log_E)).reshape(shape)
    alpha = numpy.minimum(numpy.ones_like(log_v), numpy.maximum(numpy.zeros_like(log_v), 
        numpy.log((rho / (10*min_rho)).value_in(units.none)))).reshape(shape)
    
    rgba = numpy.concatenate((red, green, blue, alpha), axis = 2)
    
    pyplot.figure(figsize = (image_size[0]/100.0, image_size[1]/100.0), dpi = 100)
    im = pyplot.figimage(rgba, origin='lower')
    
    pyplot.savefig(figname, transparent=True, dpi = 100)
    print "\nHydroplot was saved to: ", figname
    pyplot.close()


if __name__ == "__main__":
    print "Test run to mimic a supernova in SPH"
    print
    print "Details:"
    print "First a high-mass star is evolved up to the super giant phase using MESA. " \
        "Then it is converted to SPH particles with the convert_stellar_model_to_SPH " \
        "procedure (with a non-SPH 'core' particle). Finally the internal energies of " \
        "the inner particles are increased, such that the star gains the 10^51 ergs " \
        "released in supernova explosions."
    print
    run_supernova()
