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

def simulate_evolution_tracks(masses = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] | units.MSun, \
    name_of_the_figure = "HR_evolution_tracks.png", stellar_evolution_code=1):
    """
    For every mass in the `masses' array, a stellar evolution track across the Hertzsprung-Russell
    diagram will be calculated and plotted. Each star will be created, evolved and removed one by 
    one. This is only necessary because the time span of each track is different (a solar mass star
    evolution track takes billions of years, but we don't want to also evolve high mass stars for 
    billions of years) In most applications the stars have to be evolved up to a common end time,
    which can be more easily accomplished by creating an array (stars = core.Stars(number_of_stars))
    and using evolve_model(end_time = ...).
    """
    number_of_stars=len(masses)
    all_tracks_luminosity = []
    all_tracks_temperature = []
    all_tracks_stellar_type = []
    
    if stellar_evolution_code == 1:
        stellar_evolution = SSE()
    elif stellar_evolution_code == 2:
        stellar_evolution = EVtwin()
    elif stellar_evolution_code == 3:
        stellar_evolution = MESA()
#        stellar_evolution.parameters.metallicity = 0.0 | units.none
        stellar_evolution.parameters.max_iter_stop_condition = 10000 | units.none
    else:
        print "Unknown stellar_evolution_code: ", stellar_evolution_code
        return
    stellar_evolution.initialize_module_with_current_parameters() 
    
    print "The evolution across the Hertzsprung-Russell diagram of ", str(number_of_stars), \
        " stars with\nvarying masses will be simulated..."
    
    for j in range(number_of_stars):
        stars = core.Stars(1)
        star=stars[0]
        star.mass = masses[j]
        star.radius = 0.0 | units.RSun
        stars.add_calculated_attribute("temperature",calculate_effective_temperature)
        print "Created new star with mass: ", star.mass

        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        from_code_to_model = stellar_evolution.particles.new_channel_to(stars)
        from_code_to_model.copy()
    
        luminosity_at_time   = [] | units.LSun
        temperature_at_time  = [] | units.K
        stellar_type_at_time = [] | units.stellar_type
        
        stopped_evolving = False
#       Evolve this star until it changes into a compact stellar remnant (white dwarf, neutron star, or black hole)
        while star.stellar_type.value_in(units.stellar_type) < 10 and not stopped_evolving:
            luminosity_at_time.append(star.luminosity)
            temperature_at_time.append(star.temperature)
            stellar_type_at_time.append(star.stellar_type)
            previous_age = star.age
            try:
                stellar_evolution.evolve_model()
                from_code_to_model.copy()
                stopped_evolving = (star.age == previous_age) # Check whether the age has stopped increasing
            except Exception as ex:
                print str(ex)
                stopped_evolving = True
        if stopped_evolving: print "Age did not increase during timestep. Aborted evolving..."
        print " ... evolved model to t = " + str(star.age.as_quantity_in(units.Myr))
        print "Star has now become a: ", star.stellar_type, "(stellar_type: "+str(star.stellar_type.value_in(units.stellar_type))+")"
        print
        all_tracks_luminosity.append(luminosity_at_time)
        all_tracks_temperature.append(temperature_at_time)
        all_tracks_stellar_type.append(stellar_type_at_time)
        
#       Remove the star before creating the next one. See comments at the top.
        stellar_evolution.particles.remove_particles(stars)

    stellar_evolution.print_refs()
    del stellar_evolution
    
    plot_HR_diagram(number_of_stars, masses, all_tracks_luminosity, all_tracks_temperature, 
        all_tracks_stellar_type, name_of_the_figure, stellar_evolution_code)
    
    

def plot_HR_diagram(number_of_stars, masses, all_tracks_luminosity, all_tracks_temperature, 
        all_tracks_stellar_type, name_of_the_figure, stellar_evolution_code):
    if HAS_MATPLOTLIB:
        print "Plotting the data..."
        pyplot.figure(figsize = (7, 8))
        pyplot.suptitle('Hertzsprung-Russell diagram', fontsize=16)
        if stellar_evolution_code == 1:
            pyplot.title('Evolutionary tracks were simulated using the SSE package\n'
                '(Hurley J.R., Pols O.R., Tout C.A., 2000, MNRAS, 315, 543)', 
                fontsize=12)
        elif stellar_evolution_code == 2:
            pyplot.title('Evolutionary tracks were simulated using the EVtwin package', 
                fontsize=12)
        elif stellar_evolution_code == 3:
            pyplot.title('Evolutionary tracks were simulated using the MESA package', 
                fontsize=12)
        pyplot.xlabel('Effective Temperature (K)')
        pyplot.ylabel('Luminosity (solar luminosity)')

#       Define some strings for line formatting (colors, symbols, etc.), used recurrently when many stars are simulated
        plot_format_strings_lines=["r-","y-","c-","b-","m-"]
        len_fmt_str_lin=len(plot_format_strings_lines)
        plot_format_strings_symbols=["r^","y^","c^","b^","m^","rs","ys","cs","bs","ms"]
        len_fmt_str_sym=len(plot_format_strings_symbols)
        
        for j in range(number_of_stars):
#           Plot track of the current star j
            x_values = all_tracks_temperature[j].value_in(units.K)
            y_values = all_tracks_luminosity[j].value_in(units.LSun)
            pyplot.loglog(x_values, y_values, plot_format_strings_lines[j%len_fmt_str_lin])
            
#           Plot symbols whenever this star has switched to the next stellar evolution phase
            x_values = []
            y_values = []
            text_values = []
            previous_type = 15 | units.stellar_type
            for i, type in enumerate(all_tracks_stellar_type[j]):
                if not type == previous_type:
                    x_values.append(all_tracks_temperature[j][i].value_in(units.K))
                    y_values.append(all_tracks_luminosity[j][i].value_in(units.LSun))
                    text_values.append(all_tracks_stellar_type[j][i].value_in(units.stellar_type))
                    previous_type = type

            pyplot.loglog(x_values, y_values, plot_format_strings_symbols[j%len_fmt_str_sym])
            text_offset_factor_x=1.05
            text_offset_factor_y=0.6
            for i, phase in enumerate(text_values):
                pyplot.annotate(str(int(phase)), xy=(x_values[i],y_values[i]), \
                    xytext=(x_values[i]*text_offset_factor_x,y_values[i]*text_offset_factor_y))
            text_offset_factor_x=1.1
            text_offset_factor_y=0.9
            pyplot.annotate(str(masses[j]),xy=(x_values[0],y_values[0]), \
                xytext=(x_values[0]*text_offset_factor_x,y_values[0]*text_offset_factor_y), \
                color='g', horizontalalignment='right')
                
        pyplot.axis([300000., 2500., 1.e-2, 1.e6])
#       Or use these axes to also view neutron stars and black holes:
#        pyplot.axis([1.e7, 2500., 1.e-11, 1.e6])
        pyplot.savefig(name_of_the_figure)
        
        print
        print "Meaning of the stellar evolution phase markers (black numbers):"
        for i in range(16):
            print str(i)+": ", (i | units.stellar_type)
        
    print
    print "All done!"        

def calculate_effective_temperature(luminosity,radius):
    Stefan_Boltzmann_constant = 5.670400e-8 | units.J * units.s**-1 * units.m**-2 * units.K**-4
    return ((luminosity/(4*numpy.pi*Stefan_Boltzmann_constant*radius**2))**.25).in_(units.K)

def test_simulate_one_star():
    assert is_mpd_running()
    test_results_path = path_to_test_results.get_path_to_test_results()
    output_file = os.path.join(test_results_path, "HR_evolution_tracks.png")
    simulate_evolution_tracks([20.0] | units.MSun, name_of_the_figure=output_file, 
        stellar_evolution_code = 1)
    
if __name__ == '__main__':
    if not is_mpd_running():
        print "There is no mpd server running. Please do 'mpd &' first."
        sys.exit()
    if len(sys.argv) == 1:
        simulate_evolution_tracks()
    elif len(sys.argv) == 2:
        simulate_evolution_tracks(name_of_the_figure = sys.argv[1])
    else:
        if sys.argv[2] == "evtwin":
            simulate_evolution_tracks(name_of_the_figure = sys.argv[1], 
                stellar_evolution_code = 2)
        elif sys.argv[2] == "mesa":
            simulate_evolution_tracks(name_of_the_figure = sys.argv[1], 
                stellar_evolution_code = 3)
        else:
            print "Unknown option:", sys.argv[2]
            print "For EVtwin use option ' evtwin'."
            print "For MESA use option ' mesa'. Now using SSE..."
            simulate_evolution_tracks(name_of_the_figure = sys.argv[1], 
                stellar_evolution_code = 1)
