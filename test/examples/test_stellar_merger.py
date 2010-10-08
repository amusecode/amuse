import os.path
from amuse.test.amusetest import get_path_to_results
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles, ParticlesSuperset
from amuse.support.units import units, constants
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.exceptions import AmuseException
from amuse.legacy.mesa.interface import MESA
from amuse.legacy.gadget2.interface import Gadget2
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from optparse import OptionParser
import numpy

def slowtest1():
    head_on_stellar_merger()

def head_on_stellar_merger(
        masses = [1.0, 1.0] | units.MSun, 
        star_age = 10.0 | units.Gyr, 
        initial_separation = 4.0 | units.RSun, 
        initial_speed = 100.0 | units.km / units.s, 
        number_of_sph_particles = 300, 
        t_end = 2.0e4 | units.s,
        maximally_evolved_stars = False
    ):
    n_string = "n" + ("%1.0e"%(number_of_sph_particles)).replace("+0","").replace("+","")
    t_end_string = "t" + ("%1.0e"%(t_end.value_in(units.s))).replace("+0","").replace("+","")
    
    try:
        stellar_evolution = MESA()
    except:
        print "MESA was not built. Returning."
        return
    stars =  Particles(2)
    stars.mass = masses
    stellar_evolution.initialize_module_with_default_parameters() 
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
    
    composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
    outer_radii = stellar_evolution.particles[0].get_radius_profile()
    outer_radii.prepend(0.0 | units.m)
    midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
    temperature = stellar_evolution.particles[0].get_temperature_profile()
    mu          = stellar_evolution.particles[0].get_mu_profile()
    specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg)
    
    print "Creating initial conditions from a MESA stellar evolution model:"
    print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
    sph_particles_1 = convert_stellar_model_to_SPH(
        stellar_evolution.particles[0], 
        number_of_sph_particles, 
        seed=12345,
        mode = "scaling method"
    )
    print stars.mass[1], "star consisting of", number_of_sph_particles, "particles."
    sph_particles_2 = convert_stellar_model_to_SPH(
        stellar_evolution.particles[1], 
        number_of_sph_particles, 
        seed=12345,
        mode = "scaling method"
    )
    sph_particles_2.x  += initial_separation + stellar_evolution.particles.radius.sum()
    sph_particles_1.vx += initial_speed
    stellar_evolution.stop()
    all_sph_particles = ParticlesSuperset([sph_particles_1, sph_particles_2])
    
    unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, masses.sum(), t_end)
    hydro_legacy_code = Gadget2(unit_converter)
    hydro_legacy_code.gas_particles.add_particles(all_sph_particles)
    
    times = [] | units.Myr
    kinetic_energies =   [] | units.J
    potential_energies = [] | units.J
    thermal_energies =   [] | units.J
    
    print "Evolving to:", t_end
    n_steps = 100
    for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
        hydro_legacy_code.evolve_model(time)
        times.append(time)
        kinetic_energies.append(   hydro_legacy_code.kinetic_energy)
        potential_energies.append( hydro_legacy_code.potential_energy)
        thermal_energies.append(   hydro_legacy_code.thermal_energy)
    
    hydro_legacy_code.gas_particles.new_channel_to(all_sph_particles).copy_attributes(
        ['mass', 'rho', 'x','y','z', 'vx','vy','vz', 'u'])
    center_of_mass = all_sph_particles.center_of_mass().as_quantity_in(units.RSun)
    center_of_mass_velocity = all_sph_particles.center_of_mass_velocity().as_quantity_in(units.km / units.s)
    print
    print "center_of_mass:", center_of_mass
    print "center_of_mass_velocity:", center_of_mass_velocity
    all_sph_particles.position -= center_of_mass
    sph_midpoints = all_sph_particles.position.lengths()
    
    base_output_file_name = os.path.join(get_path_to_results(), "stellar_merger_"+n_string+"_"+t_end_string)
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
