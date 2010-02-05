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
from amuse.legacy.sse.muse_stellar_mpi import SSE
from amuse.legacy.evtwin.interface import EVtwin
from amuse.legacy.support.core import is_mpd_running
from support import path_to_test_results
from amuse.experiments.salpeter import SalpeterIMF

def simulate_stellar_evolution(number_of_stars = 1000, end_time = 1000.0 | units.Myr, \
    name_of_the_figure = "cluster_HR_diagram.png", use_SSE=True):
    """
    A cluster of stars will be created, with masses following a Salpeter IMF. The stellar
    evolution will be simulated using any of the legacy codes SSE, EVtwin, ...
    Finally, a Hertzsprung-Russell diagram will be produced for the final state of the 
    simulation.
    """
    print "The evolution of ", str(number_of_stars), " stars will be ", \
        "simulated until t=", str(end_time), "..."
    
    if use_SSE:
        print "Using SSE legacy code for stellar evolution."
        stellar_evolution = SSE()
    else:
        print "Using EVtwin legacy code for stellar evolution."
        stellar_evolution = EVtwin()
        if number_of_stars > stellar_evolution.parameters.maximum_number_of_stars.value_in(units.none):
            stellar_evolution.parameters.maximum_number_of_stars = (number_of_stars | units.none)
            print "You're simulating a large number of stars with EVtwin. This may be not", \
                " such a good idea..."
    stellar_evolution.initialize_module_with_current_parameters()
    
    print "Deriving a set of ", str(number_of_stars), " random masses ", \
        "following a Salpeter IMF between 0.1 and 125 MSun (alpha = -2.35)."
    initial_mass_function = SalpeterIMF()
    total_mass, salpeter_masses = initial_mass_function.next_set(number_of_stars)
    
    print "Initializing the particles"
    stars = core.Stars(number_of_stars)
    stars.mass = salpeter_masses
    stars.radius = 0.0 | units.RSun
    stars.add_calculated_attribute("temperature",calculate_effective_temperature)
    stellar_evolution.setup_particles(stars)
    stellar_evolution.initialize_stars()
    from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(stars)
    from_stellar_evolution_to_model.copy()

    print "Start evolving..."
    stellar_evolution.evolve_model(end_time)
    from_stellar_evolution_to_model.copy()

    print "Evolved model successfully."
    temperatures = stars.temperature
    luminosities = stars.luminosity
    del stellar_evolution
    
    plot_HR_diagram(temperatures, luminosities, name_of_the_figure, end_time, use_SSE)
    
    

def plot_HR_diagram(temperatures, luminosities, name_of_the_figure, end_time, use_SSE):
    if HAS_MATPLOTLIB:
        print "Plotting the data..."
        number_of_stars=len(temperatures)
        pyplot.figure(figsize = (7, 8))
        pyplot.suptitle('Hertzsprung-Russell diagram', fontsize=16)
        if use_SSE:
            pyplot.title('Stellar evolution was simulated using the SSE package\n(Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543)', \
                fontsize=12)
        else:
            pyplot.title('Stellar evolution was simulated using the EVtwin package', fontsize=12)
        pyplot.xlabel('Effective Temperature (K)')
        pyplot.ylabel('Luminosity (solar luminosity)')
        
        pyplot.loglog(temperatures.value_in(units.K), luminosities.value_in(units.LSun), "ro")


#            pyplot.loglog(x_values, y_values, plot_format_strings_symbols[j%len_fmt_str_sym])
#            text_offset_factor_x=1.05
#            text_offset_factor_y=0.6
#            for i, phase in enumerate(text_values):
#                pyplot.annotate(str(int(phase)), xy=(x_values[i],y_values[i]), \
#                    xytext=(x_values[i]*text_offset_factor_x,y_values[i]*text_offset_factor_y))
#            text_offset_factor_x=1.1
#            text_offset_factor_y=0.9
#            pyplot.annotate(str(masses[j]),xy=(x_values[0],y_values[0]), \
#                xytext=(x_values[0]*text_offset_factor_x,y_values[0]*text_offset_factor_y), \
#                color='g', horizontalalignment='right')
        
        xmin=20000.0
        xmax=2500.0
        ymin=1.e-4
        ymax=1.e4
        
        pyplot.text(xmin*.75,ymax*0.1,str(number_of_stars)+" stars\nt="+str(end_time))
        pyplot.axis([xmin, xmax, ymin, ymax])
        pyplot.savefig(name_of_the_figure)
        
    print
    print "All done!"        

def calculate_effective_temperature(luminosity,radius):
    Stefan_Boltzmann_constant = 5.670400e-8 | units.J * units.s**-1 * units.m**-2 * units.K**-4
    return ((luminosity/(4*numpy.pi*Stefan_Boltzmann_constant*radius**2))**.25).in_(units.K)

def test_simulate_one_star():
    assert is_mpd_running()
    test_results_path = path_to_test_results.get_path_to_test_results()
    output_file = os.path.join(test_results_path, "cluster_HR_diagram.png")
    simulate_stellar_evolution(100, end_time = 2.0 | units.Myr, \
        name_of_the_figure=output_file, use_SSE=True)
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        simulate_stellar_evolution()
    elif len(sys.argv) == 2:
        simulate_stellar_evolution(name_of_the_figure = sys.argv[1])
    elif len(sys.argv) == 3:
        simulate_stellar_evolution(number_of_stars=int(sys.argv[1]),end_time=(float(sys.argv[2]) | units.Myr))
    elif len(sys.argv) == 4:
        simulate_stellar_evolution(name_of_the_figure = sys.argv[1], \
            number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr))
    else:
        if sys.argv[4] == "evtwin":
            simulate_stellar_evolution(name_of_the_figure = sys.argv[1], \
                number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr), use_SSE = False)
        else:
            print "Unknown option:", sys.argv[4]
            print "For EVtwin use option ' evtwin'. Now using SSE..."
            simulate_stellar_evolution(name_of_the_figure = sys.argv[1], \
                number_of_stars=int(sys.argv[2]),end_time=(float(sys.argv[3]) | units.Myr))
