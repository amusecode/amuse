import sys
import numpy 
import os

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.units import units
from amuse.support.data import core
from amuse.legacy.sse.interface import SSE
from amuse.legacy.evtwin.interface import EVtwin
from amuse.legacy.mesa.interface import MESA
from amuse.legacy.support.core import is_mpd_running
import path_to_test_results
from amuse.ext.salpeter import SalpeterIMF

def simulate_stellar_evolution(number_of_stars = 1000, end_time = 1000.0 | units.Myr, \
    name_of_the_figure = "cluster_HR_diagram.png", stellar_evolution_code=1):
    """
    A cluster of stars will be created, with masses following a Salpeter IMF. The stellar
    evolution will be simulated using any of the legacy codes (SSE, EVtwin, or MESA).
    Finally, a Hertzsprung-Russell diagram will be produced for the final state of the 
    simulation.
    """
    print "The evolution of ", str(number_of_stars), " stars will be ", \
        "simulated until t=", str(end_time), "..."
    
    if stellar_evolution_code == 1:
        print "Using SSE legacy code for stellar evolution."
        stellar_evolution = SSE()
    elif stellar_evolution_code == 2:
        print "Using EVtwin legacy code for stellar evolution."
        stellar_evolution = EVtwin()
        if number_of_stars > stellar_evolution.parameters.maximum_number_of_stars.value_in(units.none):
            stellar_evolution.parameters.maximum_number_of_stars = (number_of_stars | units.none)
            print "You're simulating a large number of stars with EVtwin. This may be not", \
                "such a good idea..."
    elif stellar_evolution_code == 3:
        print "Using MESA legacy code for stellar evolution."
        stellar_evolution = MESA()
        if number_of_stars > 10:
            print "You're simulating a large number of stars with MESA. This may be not", \
                "such a good idea..."
    else:
        print "Unknown stellar_evolution_code: ", stellar_evolution_code
        return
    stellar_evolution.initialize_module_with_current_parameters()
    
    print "Deriving a set of ", str(number_of_stars), " random masses", \
        "following a Salpeter IMF between 0.1 and 125 MSun (alpha = -2.35)."
    initial_mass_function = SalpeterIMF()
    total_mass, salpeter_masses = initial_mass_function.next_set(number_of_stars)
    
    print "Initializing the particles"
    stars = core.Stars(number_of_stars)
    stars.mass = salpeter_masses
    print stars
    stars = stellar_evolution.particles.add_particles(stars)
    stellar_evolution.initialize_stars()

    print "Start evolving..."
    stellar_evolution.evolve_model(end_time)

    print "Evolved model successfully."
    temperatures = stars.temperature
    luminosities = stars.luminosity
    stellar_evolution.print_refs()
    stellar_evolution.stop()
    
    plot_HR_diagram(temperatures, luminosities, name_of_the_figure, end_time, 
        stellar_evolution_code)
    print "All done!"

def plot_HR_diagram(temperatures, luminosities, name_of_the_figure, end_time, 
        stellar_evolution_code):
    if HAS_MATPLOTLIB:
        print "Plotting the data..."
        number_of_stars=len(temperatures)
        pyplot.figure(figsize = (7, 8))
        pyplot.suptitle('Hertzsprung-Russell diagram', fontsize=16)
        if stellar_evolution_code == 1:
            pyplot.title('Stellar evolution was simulated using the SSE package\n'+
                '(Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543)', fontsize=12)
        elif stellar_evolution_code == 2:
            pyplot.title('Stellar evolution was simulated using the EVtwin package', fontsize=12)
        elif stellar_evolution_code == 3:
            pyplot.title('Stellar evolution was simulated using the MESA package', fontsize=12)
        pyplot.xlabel(r'T$_{\rm eff}$ (K)')
        pyplot.ylabel(r'Luminosity (L$_\odot$)')
        pyplot.loglog(temperatures.value_in(units.K), luminosities.value_in(units.LSun), "ro")

        xmin, xmax = 20000.0, 2500.0
        ymin, ymax = 1.e-4, 1.e4
        pyplot.text(xmin*.75,ymax*0.1,str(number_of_stars)+" stars\nt="+str(end_time))
        pyplot.axis([xmin, xmax, ymin, ymax])
        pyplot.savefig(name_of_the_figure)
    else:
        print "Unable to produce plot: couldn't find matplotlib."
        
def test_simulate_short():
    assert is_mpd_running()
    test_results_path = path_to_test_results.get_path_to_test_results()
    output_file = os.path.join(test_results_path, "cluster_HR_diagram.png")
    simulate_stellar_evolution(100, end_time = 2.0 | units.Myr,
        name_of_the_figure=output_file, stellar_evolution_code=1)
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        simulate_stellar_evolution()
    elif len(sys.argv) == 2:
        simulate_stellar_evolution(name_of_the_figure = sys.argv[1])
    elif len(sys.argv) == 3:
        simulate_stellar_evolution(number_of_stars=int(sys.argv[1]),
            end_time=(float(sys.argv[2]) | units.Myr))
    elif len(sys.argv) == 4:
        simulate_stellar_evolution(name_of_the_figure = sys.argv[1],
            number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr))
    else:
        if sys.argv[4] == "evtwin":
            simulate_stellar_evolution(name_of_the_figure = sys.argv[1],
                number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr), 
                stellar_evolution_code=2)
        elif sys.argv[4] == "mesa":
            simulate_stellar_evolution(name_of_the_figure = sys.argv[1],
                number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr), 
                stellar_evolution_code=3)
        else:
            print "Unknown option:", sys.argv[4]
            print "For EVtwin use option ' evtwin'."
            print "For MESA use option ' mesa'. Now using SSE..."
            simulate_stellar_evolution(name_of_the_figure = sys.argv[1], \
                number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr))
