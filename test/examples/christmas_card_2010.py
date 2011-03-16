import os.path
from amuse.test.amusetest import get_path_to_results
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles, ParticlesSuperset, Grid
from amuse.support.units import units, constants
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.exceptions import AmuseException
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
import numpy


def head_on_stellar_merger(
        masses = [0.3, 3.0] | units.MSun, 
        star_age = 310.0 | units.Myr, 
        initial_separation = 4.0 | units.RSun, 
        angle = numpy.pi / 3,
        initial_speed = 3000.0 | units.km / units.s, 
        initial_speed_perpendicular = 30.0 | units.km / units.s, 
        number_of_sph_particles = 50000, 
        t_end = 1.0e4 | units.s,
        sph_code = Fi,
    ):
    """
    masses: Mass of the two stars
    star_age: Initial age of the stars
    number_of_sph_particles: Total number of particles of both stars, divided according to their masses
    t_end: (Physical, not computational) duration of the hydrodynamics simulation
    sph_code: Code to use for the hydrodynamics simulation
    """
    
    # Convert some of the input parameters to string, for use in output file names:    
    n_string = "n" + ("%1.0e"%(number_of_sph_particles)).replace("+0","").replace("+","")
    t_end_string = "t" + ("%1.0e"%(t_end.value_in(units.s))).replace("+0","").replace("+","")
    
    base_output_file_name = os.path.join(get_path_to_results(), "stellar_merger_"+n_string+"_"+t_end_string)
    
    stars =  Particles(2)
    stars.mass = masses
    try:
        stellar_evolution = MESA()
    except:
        print "MESA was not built. Returning."
        return
    stellar_evolution.initialize_module_with_current_parameters() 
    stellar_evolution.particles.add_particles(stars)
    stellar_evolution.commit_particles()
    print "Evolving stars with MESA..."
    stellar_evolution.evolve_model(star_age)
    
    number_of_sph_particles_1 = int(round(number_of_sph_particles * 
        (stellar_evolution.particles[0].mass / stellar_evolution.particles.mass.sum()).value_in(units.none)))
    number_of_sph_particles_2 = number_of_sph_particles - number_of_sph_particles_1
    print "Creating initial conditions from a MESA stellar evolution model:"
    print stellar_evolution.particles[0].mass, "star consisting of", number_of_sph_particles_1, "particles."
    sph_particles_1 = convert_stellar_model_to_SPH(
        stellar_evolution.particles[0], 
        number_of_sph_particles_1, 
        seed=12345
    )
    print stellar_evolution.particles[1].mass, "star consisting of", number_of_sph_particles_2, "particles."
    sph_particles_2 = convert_stellar_model_to_SPH(
        stellar_evolution.particles[1], 
        number_of_sph_particles_2
    )
    
    initial_separation += stellar_evolution.particles.radius.sum()
    sph_particles_2.x  += numpy.cos(angle) * initial_separation
    sph_particles_2.y  += numpy.sin(angle) * initial_separation
    sph_particles_1.vx += numpy.cos(angle) * initial_speed - numpy.sin(angle) * initial_speed_perpendicular
    sph_particles_1.vy += numpy.cos(angle) * initial_speed_perpendicular + numpy.sin(angle) * initial_speed
    view = [-0.5, 0.5, -0.5, 0.5] * (initial_separation + stellar_evolution.particles.radius.sum())
    stellar_evolution.stop()
    
    all_sph_particles = ParticlesSuperset([sph_particles_1, sph_particles_2])
    all_sph_particles.move_to_center()
    
    unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, constants.G, t_end)
    hydro_legacy_code = sph_code(unit_converter)
    n_steps = 100
    hydro_legacy_code.parameters.n_smooth = 96 | units.none
    try:
        hydro_legacy_code.parameters.timestep = t_end / n_steps
    except Exception as exc:
        if not "parameter is read-only" in str(exc): raise
    hydro_legacy_code.gas_particles.add_particles(all_sph_particles)
    
    print "Evolving to t =", t_end, " (using", sph_code.__name__, "SPH code)."
    for time, i_step in [(i*t_end/n_steps, i) for i in range(1, n_steps+1)]:
        hydro_legacy_code.evolve_model(time)
        if not i_step % 4:
            hydro_plot(
                view,
                hydro_legacy_code,
                (300, 300),
                base_output_file_name + "_hydro_image{0:=03}.png".format(i_step)
            )
    hydro_legacy_code.stop()
    print "All done!\n"


def hydro_plot(view, hydro_code, image_size, figname):
    """
    view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    hydro_code: hydrodynamics code in which the gas to be plotted is defined
    image_size: size of the output image in pixels (x, y)
    """
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
    min_rho = 3.0e-4 | units.g / units.cm**3
    max_rho = 0.1 | units.g / units.cm**3
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
    print "Running the simulation that formed the basis of the christmas card of Leiden Observatory of 2010."
    print
    print "Details:"
    print "The ornaments are the result of a smoothed particle simulation " \
        "with 50000 equal mass particles of a 310 Myr star of 0.3 solar mass, " \
        "which is ejected from a distance of 4 solar radii (left) " \
        "with a velocity of 3000 km/s into a 3.0 solar mass star at an " \
        "age of 310 Myr. The calculation was performed using the AMUSE " \
        "(amusecode.org) software environment in which the stars were evolved " \
        "using MESA to an age of 310 Myr before the encounter was performed " \
        "using Fi. Each ornament, generated using pyplot, is a snapshot from the " \
        "simulation, from top left to bottom right. The peak is created from a " \
        "blend of all snapshots. The colors of the ornaments are, red: log " \
        "of the density, green: log of the speed and for blue we used " \
        "the log of the specific internal energy."
    print
    if HAS_MATPLOTLIB:
        head_on_stellar_merger()
    else:
        print "matplotlib is not installed. Install it in the site-packages folder of your Python installation. Returning."
