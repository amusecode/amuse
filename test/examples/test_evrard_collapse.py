import os.path
from optparse import OptionParser

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, xlabel, ylabel
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.test.amusetest import get_path_to_results
from amuse.support.units import units, constants, nbody_system
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.ext.evrard_test import new_evrard_gas_sphere
from amuse.support.codes.core import is_mpd_running
from amuse.legacy.gadget2.interface import Gadget2
from amuse.legacy.fi.interface import Fi

usage = """\
usage: amuse.sh %prog [options]

This example script will simulate the spherical collapse of 
an initially-cold adiabatic gas cloud, with a radial density 
profile proportional to 1/r.
"""

labels = [['K', 'V', 'U', 'E'], ['K', 'V', 'U', 'E']]
fi_gadget_labels = [['K_fi', 'V_fi', 'U_fi', 'E_fi'], ['K_gadget', 'V_gadget', 'U_gadget', 'E_gadget']]

convert_generic_units = ConvertBetweenGenericAndSiUnits(1.0 | units.kpc, 1.0e10 | units.MSun, constants.G)
convert_nbody_units   = nbody_system.nbody_to_si(       1.0 | units.kpc, 1.0e10 | units.MSun)

def run_evrard(
        hydro_legacy_codes, 
        number_of_particles, 
        random_seed = None, 
        name_of_the_figure = "evrard_collapse_test.png"):
    
    print "\nThe spherical collapse of an initially-cold adiabatic gas cloud,\n", \
        "consisting of ", str(number_of_particles), "particles will be simulated...\n"
    
    t_end = 3.0 * convert_nbody_units.to_si(1.0 | nbody_system.time) # (3 natural timescales)
    print "Evolving to (3 natural timescales): ", t_end.as_quantity_in(units.Myr)
    n_steps = 100
    
    gas = new_evrard_gas_sphere(number_of_particles, convert_nbody_units, do_scale = True, seed = random_seed)
    gas.h_smooth = 0.01 | units.kpc
    
    times = [] | units.Myr
    kinetic_energies =   []
    potential_energies = []
    thermal_energies =   []
    
    for hydro_legacy_code in hydro_legacy_codes:
        try:
            hydro_legacy_code.parameters.timestep = t_end / n_steps
        except Exception as exc:
            if not "parameter is read-only" in str(exc): raise
        
        hydro_legacy_code.gas_particles.add_particles(gas)
        
        kinetic_energies.append([] | units.J)
        potential_energies.append([] | units.J)
        thermal_energies.append([] | units.J)
    
    for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
        for i, hydro_legacy_code in enumerate(hydro_legacy_codes):
            hydro_legacy_code.evolve_model(time)
            kinetic_energies[i].append(   hydro_legacy_code.kinetic_energy)
            potential_energies[i].append( hydro_legacy_code.potential_energy)
            thermal_energies[i].append(   hydro_legacy_code.thermal_energy)
        times.append(time)
    hydro_legacy_code.stop()
    energy_plot(times, kinetic_energies, potential_energies, thermal_energies, name_of_the_figure)
    print "All done!\n"

def energy_plot(time, E_kin_list, E_pot_list, E_therm_list, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    for i, (E_kin, E_pot, E_therm) in enumerate(zip(E_kin_list, E_pot_list, E_therm_list)):
        plot(time, E_kin.as_quantity_in(units.erg), label=labels[i][0])
        plot(time, E_pot, label=labels[i][1])
        plot(time, E_therm, label=labels[i][2])
        plot(time, E_kin+E_pot+E_therm, label=labels[i][3])
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(prop={'size':"x-small"}, loc=3)
    pyplot.savefig(figname)
    print "\nPlot of energy evolution was saved to: ", figname
    pyplot.close()

class InstantiateCode(object):
    
    def fi(self):
        return Fi(convert_nbody_units)
    
    def gadget(self):
        return Gadget2(convert_generic_units)

    def new_code(self, name_of_the_code):
        if hasattr(self, name_of_the_code):
            return getattr(self, name_of_the_code)()
        else:
            raise Exception("Cannot instantiate code with name '{0}'".format(name_of_the_code))

def new_code(name_of_the_code):
    return InstantiateCode().new_code(name_of_the_code)

def test_evrard_fi_short():
    assert is_mpd_running()
    run_evrard(
        [new_code("fi")],
        100,
        random_seed = 12345, 
        name_of_the_figure = os.path.join(get_path_to_results(), "evrard_test_fi_100.png")
    )

def test_evrard_gadget_short():
    assert is_mpd_running()
    run_evrard(
        [new_code("gadget")],
        100,
        random_seed = 12345, 
        name_of_the_figure = os.path.join(get_path_to_results(), "evrard_test_gadget_100.png")
    )

def test_evrard_compare_short():
    assert is_mpd_running()
    global labels
    labels = fi_gadget_labels
    run_evrard(
        [new_code("fi"), new_code("gadget")],
        100,
        random_seed = 12345, 
        name_of_the_figure = os.path.join(get_path_to_results(), "evrard_test_compare_100.png")
    )

def slowtest_evrard_fi():
    assert is_mpd_running()
    run_evrard(
        [new_code("fi")],
        10000,
        random_seed = 12345, 
        name_of_the_figure = os.path.join(get_path_to_results(), "evrard_test_fi_10000.png")
    )

def slowtest_evrard_gadget():
    assert is_mpd_running()
    run_evrard(
        [new_code("gadget")],
        10000,
        random_seed = 12345, 
        name_of_the_figure = os.path.join(get_path_to_results(), "evrard_test_gadget_10000.png")
    )

def slowtest_evrard_compare():
    assert is_mpd_running()
    global labels
    labels = fi_gadget_labels
    run_evrard(
        [new_code("fi"), new_code("gadget")],
        10000,
        random_seed = 12345, 
        name_of_the_figure = os.path.join(get_path_to_results(), "evrard_test_compare_10000.png")
    )

def new_commandline_option_parser():
    result = OptionParser(usage)
    result.add_option(
        "-c",
        "--code",
        choices=["fi", "gadget"],
        default="fi",
        dest="code",
        metavar="CODE",
        help="CODE to use for hydrodynamics"
    )
    result.add_option(
        "-a",
        "--all_codes",
        default=False,
        action="store_true",
        dest="all_codes",
        metavar="ALLCODE",
        help="ALLCODES flag to use all hydrodynamics solvers"
    )
    result.add_option(
        "-n",
        "--number_of_particles",
        type="int",
        default=100,
        dest="number_of_particles",
        help="Number of particles in the gas cloud"
    )
    result.add_option(
        "-S",
        "--seed",
        type="int",
        default=None,
        dest="salpeterSeed",
        help="Random number generator seed"
    )
    result.add_option(
        "-p",
        "--plot_file",
        type="string",
        default="evrard_collapse_test.png",
        dest="plot_filename",
        help="Name of the file to plot to"
    )
    return result

if __name__ == '__main__':
    if not is_mpd_running():
        print "There is no mpd server running. Please do 'mpd &' first."
        sys.exit()
    parser = new_commandline_option_parser()
    (options, arguments) = parser.parse_args()
    if arguments:
        parser.error("unknown arguments '{0}'".format(arguments))
    
    if options.all_codes:
        code = [new_code("fi"), new_code("gadget")]
        labels = fi_gadget_labels
    else:
        code = [new_code(options.code)]
    run_evrard(
        code,
        options.number_of_particles,
        random_seed = options.salpeterSeed, 
        name_of_the_figure = options.plot_filename
    )
