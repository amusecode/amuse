import numpy
import os.path
from amuse.test.amusetest import get_path_to_results
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.units import units
from amuse.units import generic_unit_system
from amuse.units import nbody_system
from amuse.units import constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.exceptions import AmuseException
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.ext.star_to_sph import *

from amuse.datamodel import Particles
from amuse.datamodel import Grid

def setup_stellar_evolution_model():
    out_pickle_file = os.path.join(get_path_to_results(),
                                   "super_giant_stellar_structure.pkl")
    if os.path.exists(out_pickle_file):
        return out_pickle_file

    print "Creating initial conditions from a MESA stellar evolution model..."

    stellar_evolution = MESA(redirection = "none")
    stars =  Particles(1)
    stars.mass = 10.0 | units.MSun
    stellar_evolution.particles.add_particles(stars)

    while stellar_evolution.particles[0].stellar_type <= 12|units.stellar_type:
        stellar_evolution.evolve_model()
    
    pickle_stellar_model(stellar_evolution.particles[0], out_pickle_file)
    stellar_evolution.stop()
    return out_pickle_file

def inject_supernova_energy(gas_particles):
    Rinner = 10|units.RSun
    inner = gas_particles.select(
        lambda pos : pos.length_squared() < Rinner**2, ["position"])
    print "Adding", (1.0e51 | units.erg) / inner.total_mass(),
    print "to each of", len(inner), "innermost particles"
    print "    of the exploding star"
    inner.u += (1.0e51 | units.erg) / inner.total_mass()

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def hydro_plot(view, hydro_code, image_size, time, figname):
    """
    Produce a series of images suitable for conversion into a movie.
    view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    hydro_code: hydrodynamics code in which the gas to be plotted is defined
    image_size: size of the output image in pixels (x, y)
    time: current hydro code time
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
    rho, rhovx, rhovy, rhovz, rhoe \
        = hydro_code.get_hydro_state_at_point(
            grid.x.reshape(size), grid.y.reshape(size), grid.z.reshape(size),
            speed, speed, speed)
    
    min_v =  800.0 | units.km / units.s
    max_v = 3000.0 | units.km / units.s
    min_rho = 3.0e-9 | units.g / units.cm**3
    max_rho = 1.0e-5 | units.g / units.cm**3
    min_E = 1.0e11 | units.J / units.kg
    max_E = 1.0e13 | units.J / units.kg
    
    v_sqr = (rhovx**2 + rhovy**2 + rhovz**2) / rho**2
    E = rhoe / rho
    log_v = numpy.log((v_sqr/min_v**2)) / numpy.log((max_v**2/min_v**2))
    log_rho = numpy.log((rho/min_rho)) / numpy.log((max_rho/min_rho))
    log_E = numpy.log((E/min_E)) / numpy.log((max_E/min_E))
    
    red   = numpy.minimum(numpy.ones_like(rho.number),
                          numpy.maximum(numpy.zeros_like(rho.number),
                                        log_rho)).reshape(shape)
    green = numpy.minimum(numpy.ones_like(rho.number),
                          numpy.maximum(numpy.zeros_like(rho.number),
                                        log_v)).reshape(shape)
    blue  = numpy.minimum(numpy.ones_like(rho.number),
                          numpy.maximum(numpy.zeros_like(rho.number),
                                        log_E)).reshape(shape)
    alpha = numpy.minimum(numpy.ones_like(log_v),
                          numpy.maximum(numpy.zeros_like(log_v), 
                                        numpy.log((rho
                                              / (10*min_rho))))).reshape(shape)
    
    rgba = numpy.concatenate((red, green, blue, alpha), axis = 2)
    
    pyplot.figure(figsize = (image_size[0]/100.0,
                             image_size[1]/100.0), dpi = 100)
    im = pyplot.figimage(rgba, origin='lower')
    
    pyplot.savefig(figname, transparent=True, dpi = 100,
                   facecolor='k', edgecolor='k')
    print "Saved hydroplot at time", time, "in file"
    print '   ', figname
    pyplot.close()

def energy_plot(time, E_kin, E_pot, E_therm, figname):
    if not HAS_MATPLOTLIB:
        return
    x_label = 'Time [hour]'
    y_label = 'Energy [foe]'
    single_frame(x_label, y_label, logx=False, logy=False,
                 xsize=14, ysize=10, ymin=-1, ymax=-1)
    cols = get_distinct(4)

    FOE = 1.e+51 | units.erg
    hour = 1|units.hour
    pyplot.plot(time/hour, E_kin/FOE, label='E_kin', c=cols[0])
    pyplot.plot(time/hour, E_pot/FOE, label='E_pot', c=cols[1])
    pyplot.plot(time/hour, E_therm/FOE, label='E_therm', c=cols[2])
    pyplot.plot(time/hour, (E_kin+E_pot+E_therm)/FOE, label='E_total',
                c=cols[3])
    pyplot.legend(loc='best')
    pyplot.savefig(figname)
    
    print '\nSaved energy evolution figure in file', figname, '\n'
    pyplot.show()
    pyplot.close()

def run_supernova():
    use_hydro_code = Gadget2
    hydro_code_options = dict(number_of_workers=3)
    number_of_sph_particles = 3000
    t_end = 1.0e4 | units.s
    
    pickle_file = setup_stellar_evolution_model()
    model = convert_stellar_model_to_SPH(None, 
                                         number_of_sph_particles, 
                                         seed = 12345,
                                         pickle_file = pickle_file,
                                         with_core_particle = True,
                                         target_core_mass = 1.4|units.MSun)
    print "model=", model.core_particle
    core, gas_without_core, core_radius \
        = model.core_particle, model.gas_particles, model.core_radius
    
    inject_supernova_energy(gas_without_core)
    
    print "\nEvolving (SPH) to:", t_end
    n_steps = 100
    
    unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun,
                                                     constants.G, t_end)
    hydro_code = use_hydro_code(unit_converter, **hydro_code_options)
    
    try:
        hydro_code.parameters.timestep = t_end / n_steps
    except Exception as exc:
        if not "parameter is read-only" in str(exc): raise
    
    hydro_code.parameters.epsilon_squared = core_radius**2
    hydro_code.parameters.n_smooth_tol = 0.01
    hydro_code.gas_particles.add_particles(gas_without_core)
    hydro_code.dm_particles.add_particle(core)
    
    times = [] | units.s
    potential_energies = [] | units.J
    kinetic_energies =   [] | units.J
    thermal_energies =   [] | units.J
    for time, i_step in [(i*t_end/n_steps, i) for i in range(0, n_steps+1)]:
        hydro_code.evolve_model(time)
        times.append(time)
        potential_energies.append(hydro_code.potential_energy)
        kinetic_energies.append(  hydro_code.kinetic_energy)
        thermal_energies.append(  hydro_code.thermal_energy)
        hydro_plot([-1.0, 1.0, -1.0, 1.0] * (350 | units.RSun),
                   hydro_code,
                   (100, 100),
                   time,
                   os.path.join(get_path_to_results(),
                            "supernova_hydro_image{0:=03}.png".format(i_step))
        )
    
    energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
                "supernova_energy_evolution.pdf")

    hydro_code.stop()
    
if __name__ == "__main__":
    print "Test run to mimic a supernova in SPH"
    print "Details:"
    print "    A high-mass star is evolved to the supergiant phase using MESA."
    print "    Then it is converted to SPH particles using", \
          "convert_stellar_model_to_SPH"
    print "    (with a non-SPH 'core' particle).", \
          "Finally the internal energies of"
    print "    the innermost particles are increased so that the star gains the"
    print "    10^51 erg released in a typical supernova explosion."

    run_supernova()
