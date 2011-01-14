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
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH, pickle_stellar_model, StellarModel2SPH
from optparse import OptionParser
import numpy

def slowtest1():
    head_on_stellar_merger()

def head_on_stellar_merger(
        masses = [0.3, 3.0] | units.MSun, 
        star_age = 310.0 | units.Myr, 
        maximally_evolved_stars = False,
        initial_separation = 4.0 | units.RSun, 
        angle = numpy.pi / 3,
        initial_speed = 3000.0 | units.km / units.s, 
        initial_speed_perpendicular = 30.0 | units.km / units.s, 
        number_of_sph_particles = 1000, 
        t_end = 1.0e4 | units.s,
        sph_code = Fi,
        steps_per_snapshot = 4,
        snapshot_size = 100,
        use_stored_stellar_models = True
    ):
    """
    masses: Mass of the two stars
    star_age: Initial age of the stars (if maximally_evolved_stars is False)
    maximally_evolved_stars: Evolve stars as far as the Stellar Evolution code can get
    number_of_sph_particles: Total number of particles of both stars, divided according to their masses
    t_end: (Physical, not computational) duration of the hydrodynamics simulation
    sph_code: Code to use for the hydrodynamics simulation
    steps_per_snapshot: A hydroplot snapshot is generated each time after this many steps (0 or None means no snapshots)
    snapshot_size: Size of the snapshot in pixels along one dimension
    use_stored_stellar_models: Flag to use previously stored stellar model files (for speed-up).
    """
    
    # Convert some of the input parameters to string, for use in output file names:    
    n_string = "n" + ("%1.0e"%(number_of_sph_particles)).replace("+0","").replace("+","")
    t_end_string = "t" + ("%1.0e"%(t_end.value_in(units.s))).replace("+0","").replace("+","")
    masses_string = ("m1_" + ("%0.3e"%(masses[0].value_in(units.MSun))).replace("+0","").replace("+","") +
        "_m2_" + ("%0.3e"%(masses[1].value_in(units.MSun))).replace("+0","").replace("+",""))
    if maximally_evolved_stars:
        star_age_string = "a_max"
    else:
        star_age_string = "a" + ("%0.3e"%(star_age.value_in(units.Myr))).replace("+0","").replace("+","")
    
    base_output_file_name = os.path.join(get_path_to_results(), "stellar_merger_"+n_string+"_"+t_end_string)
    pickle_file_1 = os.path.join(get_path_to_results(), "stellar_merger_"+masses_string+"_"+star_age_string+"_1.pkl")
    pickle_file_2 = os.path.join(get_path_to_results(), "stellar_merger_"+masses_string+"_"+star_age_string+"_2.pkl")
    
    if not use_stored_stellar_models or not (os.path.exists(pickle_file_1) and os.path.exists(pickle_file_2)):
        stars =  Particles(2)
        stars.mass = masses
        try:
            stellar_evolution = MESA()
        except:
            print "MESA was not built. Returning."
            return
        stellar_evolution.initialize_module_with_current_parameters() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.initialize_stars()
        
        if maximally_evolved_stars:
            try:
                while True:
                    stellar_evolution.evolve_model()
            except AmuseException as exception:
                print exception
        else:
            stellar_evolution.evolve_model(star_age)
        
        if os.path.exists(pickle_file_1):
            print "Could not save stellar model 1: file already exists."
        else:
            pickle_stellar_model(stellar_evolution.particles[0], pickle_file_1)
            print "Stellar model 1 saved at:", pickle_file_1
        if os.path.exists(pickle_file_2):
            print "Could not save stellar model 2: file already exists."
        else:
            pickle_stellar_model(stellar_evolution.particles[1], pickle_file_2)
            print "Stellar model 2 saved at:", pickle_file_2
        
        stellar_evolution.stop()
    
    model_1 = StellarModel2SPH(None, None, pickle_file = pickle_file_1)
    model_2 = StellarModel2SPH(None, None, pickle_file = pickle_file_2)
    model_1.unpickle_stellar_structure()
    model_2.unpickle_stellar_structure()
    composition = model_2.composition_profile
    midpoints = model_2.midpoints_profile[1:-1]
    specific_internal_energy = model_2.specific_internal_energy_profile
    
    number_of_sph_particles_1 = int(round(number_of_sph_particles * 
        (model_1.mass / (model_1.mass + model_2.mass)).value_in(units.none)))
    number_of_sph_particles_2 = number_of_sph_particles - number_of_sph_particles_1
    print "Creating initial conditions from a MESA stellar evolution model:"
    print model_1.mass, "star consisting of", number_of_sph_particles_1, "particles."
    sph_particles_1 = convert_stellar_model_to_SPH(
        None, 
        number_of_sph_particles_1, 
        seed=12345,
        mode = "scaling method",
        pickle_file = pickle_file_1
    )
    print model_2.mass, "star consisting of", number_of_sph_particles_2, "particles."
    sph_particles_2 = convert_stellar_model_to_SPH(
        None, 
        number_of_sph_particles_2, 
        mode = "scaling method",
        pickle_file = pickle_file_2
    )
    initial_separation += model_1.radius + model_2.radius
    sph_particles_2.x  += numpy.cos(angle) * initial_separation
    sph_particles_2.y  += numpy.sin(angle) * initial_separation
    sph_particles_1.vx += numpy.cos(angle) * initial_speed - numpy.sin(angle) * initial_speed_perpendicular
    sph_particles_1.vy += numpy.cos(angle) * initial_speed_perpendicular + numpy.sin(angle) * initial_speed
    view = [-0.5, 0.5, -0.5, 0.5] * (initial_separation + model_1.radius + model_2.radius)
    
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
    
    times = [] | units.Myr
    kinetic_energies =   [] | units.J
    potential_energies = [] | units.J
    thermal_energies =   [] | units.J
    
    print "Evolving to:", t_end
    for time, i_step in [(i*t_end/n_steps, i) for i in range(1, n_steps+1)]:
        hydro_legacy_code.evolve_model(time)
        times.append(time)
        kinetic_energies.append(   hydro_legacy_code.kinetic_energy)
        potential_energies.append( hydro_legacy_code.potential_energy)
        thermal_energies.append(   hydro_legacy_code.thermal_energy)
        if steps_per_snapshot and (not i_step % steps_per_snapshot):
            hydro_plot(
                view,
                hydro_legacy_code,
                (snapshot_size, snapshot_size),
                base_output_file_name + "_hydro_image{0:=03}.png".format(i_step)
            )

    
    hydro_legacy_code.gas_particles.new_channel_to(all_sph_particles).copy_attributes(
        ['mass', 'x','y','z', 'vx','vy','vz', 'u'])
    center_of_mass = all_sph_particles.center_of_mass().as_quantity_in(units.RSun)
    center_of_mass_velocity = all_sph_particles.center_of_mass_velocity().as_quantity_in(units.km / units.s)
    print
    print "center_of_mass:", center_of_mass
    print "center_of_mass_velocity:", center_of_mass_velocity
    all_sph_particles.position -= center_of_mass
    sph_midpoints = all_sph_particles.position.lengths()
    
    energy_plot(
        times, 
        kinetic_energies, potential_energies, thermal_energies, 
        base_output_file_name+"_energy_evolution.png"
    )
    thermal_energy_plot(
        times, 
        thermal_energies, 
        base_output_file_name+"_thermal_energy_evolution.png"
    )
    composition_comparison_plot(
        midpoints, composition[0], 
        sph_midpoints, all_sph_particles.h1, 
        base_output_file_name+"_composition_h1.png"
    )
    internal_energy_comparison_plot(
        midpoints, specific_internal_energy, 
        sph_midpoints, all_sph_particles.u, 
        base_output_file_name+"_new_u.png"
    )
    hydro_plot(
        [-2.0, 2.0, -2.0, 2.0] | units.RSun,
        hydro_legacy_code,
        (100, 100),
        base_output_file_name + "_hydro_image.png"
    )
    hydro_legacy_code.stop()
    print "All done!\n"


def composition_comparison_plot(radii_SE, comp_SE, radii_SPH, comp_SPH, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (7, 5))
    plot(radii_SE.as_quantity_in(units.RSun), comp_SE.value_in(units.none), 
        label='stellar evolution model')
    plot(radii_SPH, comp_SPH.value_in(units.none), 'go', label='SPH model')
    xlabel('radius')
    ylabel('mass fraction')
    pyplot.legend()
    pyplot.savefig(figname)
    print "\nPlot of composition profiles was saved to: ", figname
    pyplot.close()

def internal_energy_comparison_plot(radii_SE, u_SE, radii_SPH, u_SPH, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (7, 5))
    semilogy(radii_SE.as_quantity_in(units.RSun), u_SE, 
        label='stellar evolution model')
    semilogy(radii_SPH, u_SPH, 'go', label='SPH model')
    xlabel('radius')
    ylabel('internal energy')
    pyplot.legend()
    pyplot.savefig(figname)
    print "\nPlot of internal energy profiles was saved to: ", figname
    pyplot.close()

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

def thermal_energy_plot(time, E_therm, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    plot(time, E_therm.as_quantity_in(units.erg), label='E_therm')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(loc=3)
    pyplot.savefig(figname)
    print "\nPlot of thermal energy evolution was saved to: ", figname
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
    
    
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n",
        "--n_sph", 
        dest="number_of_sph_particles",
        type="int",
        default = 300,
        help="Number of sph particles per star in solar mass"
    )
    result.add_option(
        "-a",
        "--mass1", 
        dest="mass1",
        type="float",
        default = 1.0,
        help="Mass of the first star"
    )
    result.add_option(
        "-b",
        "--mass2", 
        dest="mass2",
        type="float",
        default = 1.0,
        help="Mass of the second star in solar mass"
    )
    result.add_option(
        "-d",
        "--distance",
        dest="initial_separation",
        type="float",
        default=4.0,
        help="Initial_separation between the two stars in solar radii"
    )
    result.add_option(
        "-v",
        "--speed",
        dest="initial_speed",
        type="float",
        default=100.0,
        help="Initial relative speed between the two stars in km/s"
    )
    result.add_option(
        "-t",
        "--end_time",
        dest="end_time",
        type="float",
        default=2.0e4,
        help="Time to evolve to in s"
    )
    result.add_option(
        "-s",
        "--t_star",
        dest="t_star",
        type="float",
        default=1.0e4,
        help="Age of both stars in Myr"
    )
    result.add_option(
        "-e",
        "--max_evolved",
        default=False,
        action="store_true",
        dest="maximally_evolved_stars",
        help="Flag to use maximally evolved stars, i.e. until MESA exits"
    )
    return result


if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    print "Simulating a head on collision between two stars."
    print "Options used:", options
    masses = [options.mass1, options.mass2] | units.MSun
    head_on_stellar_merger(
        masses = masses,
        star_age = options.t_star | units.Myr,
        initial_separation = options.initial_separation | units.RSun,
        initial_speed = options.initial_speed | units.km / units.s,
        number_of_sph_particles = options.number_of_sph_particles,
        t_end = options.end_time | units.s,
        maximally_evolved_stars = options.maximally_evolved_stars
    )
