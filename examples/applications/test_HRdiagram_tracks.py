import sys
import numpy 
import os
import warnings

from optparse import OptionParser

from amuse.units import units
from amuse.support.data import core

from amuse.community.sse.interface import SSE
from amuse.community.evtwin.interface import EVtwin
from amuse.community.evtwin2sse.interface import EVtwin2SSE
from amuse.community.mesa.interface import MESA

from amuse.community.cachedse.interface import CachedStellarEvolution

from amuse.support.codes.core import is_mpd_running
from amuse.test.amusetest import get_path_to_results

usage = """\
usage: %prog [options]
    
This script will generate HR diagram tracks for 
stars with specified a masses.
"""
def stellar_remnant_state(star):
    return 10 <= star.stellar_type.value_in(units.stellar_type) and \
        star.stellar_type.value_in(units.stellar_type) < 16

def simulate_evolution_tracks(
    stellar_evolution = SSE(),
    masses = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] | units.MSun,
    name_of_the_figure = "HR_evolution_tracks.png"
    ):
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
    stellar_evolution.commit_parameters() 
    
    print "The evolution across the Hertzsprung-Russell diagram of ", str(number_of_stars), \
        " stars with\nvarying masses will be simulated..."
    
    for j in range(number_of_stars):
        star = core.Particle()
        star.mass = masses[j]
        print "Created new star with mass: ", star.mass

        star = stellar_evolution.particles.add_particle(star)
        stellar_evolution.commit_particles()
    
        luminosity_at_time     = [] | units.LSun
        temperature_at_time     = [] | units.K
        stellar_type_at_time = [] | units.stellar_type
        
        stopped_evolving = False
#        Evolve this star until it changes into a compact stellar remnant (white dwarf, neutron star, or black hole)
        while not stellar_remnant_state(star) and not stopped_evolving:
            luminosity_at_time.append(star.luminosity)
            temperature_at_time.append(star.temperature)
            stellar_type_at_time.append(star.stellar_type)
            previous_age = star.age
            try:
                stellar_evolution.evolve_model()
                stopped_evolving = (star.age == previous_age) # Check whether the age has stopped increasing
            except Exception as ex:
                print str(ex)
                stopped_evolving = True
        if stopped_evolving:
            print "Age did not increase during timestep. Aborted evolving..."
        else:
            stellar_type_at_time.append(star.stellar_type)
            # Fudged: final stellar type annotation at previous (Teff, L);
            # BHs and neutron stars would otherwise fall off the chart.
            luminosity_at_time.append(luminosity_at_time[-1])
            temperature_at_time.append(temperature_at_time[-1])
        print " ... evolved model to t = " + str(star.age.as_quantity_in(units.Myr))
        print "Star has now become a: ", star.stellar_type, "(stellar_type: "+str(star.stellar_type.value_in(units.stellar_type))+")"
        print
        all_tracks_luminosity.append(luminosity_at_time)
        all_tracks_temperature.append(temperature_at_time)
        all_tracks_stellar_type.append(stellar_type_at_time)
        
#        Remove the star before creating the next one. See comments at the top.
        stellar_evolution.particles.remove_particle(star)

    stellar_evolution.stop()
    
    plot_HR_diagram(masses, all_tracks_luminosity, all_tracks_temperature, 
        all_tracks_stellar_type, name_of_the_figure)
    
    print "All done!"         

def plot_HR_diagram(masses, luminosity_tracks, temperature_tracks, stellar_type_tracks, plotfile):
    try:
#        This removes the need for ssh -X to be able to do plotting
        import matplotlib
        matplotlib.use("Agg", warn=False) 
        
        from matplotlib import pyplot
        print "Plotting the data..."
        pyplot.figure(figsize = (7, 8))
        pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
        pyplot.xlabel('Effective Temperature (K)')
        pyplot.ylabel('Luminosity (solar luminosity)')

#        Define some strings for line formatting (colors, symbols, etc.), used recurrently when many stars are simulated
        plot_format_strings_lines=["r-","y-","c-","b-","m-"]
        len_fmt_str_lin=len(plot_format_strings_lines)
        plot_format_strings_symbols=["r^","y^","c^","b^","m^","rs","ys","cs","bs","ms"]
        len_fmt_str_sym=len(plot_format_strings_symbols)
        
        number_of_stars = len(masses)
        for j in range(number_of_stars):
#            Plot track of the current star j
            x_values = temperature_tracks[j].value_in(units.K)
            y_values = luminosity_tracks[j].value_in(units.LSun)
            pyplot.loglog(x_values, y_values, plot_format_strings_lines[j%len_fmt_str_lin])
            
#            Plot symbols whenever this star has switched to the next stellar evolution phase
            x_values = []
            y_values = []
            text_values = []
            previous_type = 15 | units.stellar_type
            for i, type in enumerate(stellar_type_tracks[j]):
                if not type == previous_type:
                    x_values.append(temperature_tracks[j][i].value_in(units.K))
                    y_values.append(luminosity_tracks[j][i].value_in(units.LSun))
                    text_values.append(stellar_type_tracks[j][i].value_in(units.stellar_type))
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
#        Or use these axes to also view neutron stars and black holes:
#         pyplot.axis([1.e7, 2500., 1.e-11, 1.e6])
        pyplot.savefig(plotfile)
        print "Meaning of the stellar evolution phase markers (black numbers):"
        for i in range(16):
            print str(i)+": ", (i | units.stellar_type)           
    except ImportError:
        print "Unable to produce plot: couldn't find matplotlib."

class InstantiateCode(object):
    
    def sse(self, number_of_stars):
        return SSE()
    
    def evtwin(self, number_of_stars):
        result = EVtwin()
        result.initialize_code()
        if number_of_stars > result.parameters.maximum_number_of_stars:
            result.parameters.maximum_number_of_stars = number_of_stars
            warnings.warn("You're simulating a large number of stars with EVtwin. This may not be such a good idea...")
        return result
        
    def mesa(self, number_of_stars):
        result = MESA()
        result.initialize_code()
        if number_of_stars > (10 | units.none):
            warnings.warn("You're simulating a large number of stars with MESA. This may not be such a good idea...")
        if number_of_stars > (1000| units.none):
            raise Exception("You want to simulate with more than 1000 stars using MESA, this is not supported")
        return result

    def evtwin2sse(self, number_of_stars):
        result = EVtwin2SSE()
        result.initialize_code()
        # TODO add maximum_number_of_stars parameter to Evtwin2SSE
        #if number_of_stars > result.parameters.maximum_number_of_stars:
        #     result.parameters.maximum_number_of_stars = number_of_stars
        #     warnings.warn("You're simulating a large number of stars with EVtwin. This may not be such a good idea...")
        return result

    def new_code(self, name_of_the_code, number_of_stars):
        if hasattr(self, name_of_the_code):
            return getattr(self, name_of_the_code)(number_of_stars)
        else:
            raise Exception("Cannot instantiate code with name '{0}'".format(name_of_the_code))

def new_code(name_of_the_code, number_of_stars):
    return InstantiateCode().new_code(name_of_the_code, number_of_stars)

def test_simulate_one_star():
    assert is_mpd_running()
    code = new_code("sse", 1 | units.none)
    test_results_path = get_path_to_results()
    output_file = os.path.join(test_results_path, "HR_evolution_tracks.png")
    simulate_evolution_tracks(
        code,
        [20.0] | units.MSun,
        name_of_the_figure=output_file, 
    )

def new_commandline_option_parser():
    result = OptionParser(usage)
    result.add_option(
        "-c",
        "--code",
        choices=["sse", "evtwin", "mesa", "evtwin2sse"],
        default="sse",
        dest="code",
        metavar="CODE",
        help="CODE to use for stellar evolution"
    )
    result.add_option(
        "-C",
        "--cache",
        type="string",
        default=None,
        dest="cacheDir",
        help="Use/write cache from directory"
    )     
    result.add_option(
        "-p",
        "--plot_file",
        type="string",
        default="HR_evolution_tracks.png",
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
    mass_list = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] | units.MSun
    
    code = new_code(options.code, len(mass_list) | units.none)
    if not (options.cacheDir is None):
        print "Using cache directory: %s" % (options.cacheDir)
# As a special case, we use caching of the underlying models instead of the model output for EVtwin2SSE
        if (options.code == "evtwin2sse"):
            code.cache_underlying_models(options.cacheDir)
        else:
            code = CachedStellarEvolution(code, options.cacheDir)     
    simulate_evolution_tracks(
        code,
        masses = mass_list,
        name_of_the_figure=options.plot_filename
    )
