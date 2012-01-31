import sys
import numpy 
import os
import warnings

from optparse import OptionParser

from amuse.units import units
from amuse.community.sse.interface import SSE
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA
from amuse.community.cachedse.interface import CachedStellarEvolution

from amuse.test.amusetest import get_path_to_results

import numpy

from amuse import datamodel
from amuse.rfi.core import is_mpd_running
from amuse.ic.salpeter import new_salpeter_mass_distribution
usage = """\
usage: %prog [options]
	
This script will generate HR diagram for an 
evolved cluster of stars with a Salpeter mass 
distribution.
"""


def simulate_stellar_evolution(
	stellar_evolution = SSE(),
	number_of_stars = 1000, 
	end_time = 1000.0 | units.Myr, 
	name_of_the_figure = "cluster_HR_diagram.png", 
	):
    """
    A cluster of stars will be created, with masses following a Salpeter IMF. The stellar
    evolution will be simulated using any of the legacy codes (SSE, EVtwin, or MESA).
    Finally, a Hertzsprung-Russell diagram will be produced for the final state of the 
    simulation.
    """
    print "The evolution of", str(number_of_stars), "stars will be ",  \
            "simulated until t =", str(end_time), "..."
    
    stellar_evolution.commit_parameters()
    
    print ("Deriving a set of", str(number_of_stars), "random masses",  
            "following a Salpeter IMF between 0.1 and 125 MSun (alpha = -2.35).")
    
    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)
    
    print "Initializing the particles"
    stars = datamodel.Particles(number_of_stars)
    stars.mass = salpeter_masses
    print "Stars to evolve:"
    print stars
    stars = stellar_evolution.particles.add_particles(stars)
    stellar_evolution.commit_particles()

    #print stars.temperature

    print "Start evolving..."
    stellar_evolution.evolve_model(end_time)

    print "Evolved model successfully."
    temperatures = stars.temperature
    luminosities = stars.luminosity
    
    stellar_evolution.stop()
    
    plot_HR_diagram(temperatures, luminosities, name_of_the_figure, end_time)
            
    print "All done!"

def plot_HR_diagram(temperatures, luminosities, name_of_the_figure, end_time):
    try:
#		This removes the need for ssh -X to be able to do plotting
        import matplotlib
        matplotlib.use("Agg", warn=False) 
        from matplotlib import pyplot
        print "Plotting the data..."
        number_of_stars=len(temperatures)
        pyplot.figure(figsize = (7, 8))
        pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
        pyplot.xlabel(r'T$_{\rm eff}$ (K)')
        pyplot.ylabel(r'Luminosity (L$_\odot$)')
        pyplot.loglog(temperatures.value_in(units.K), luminosities.value_in(units.LSun), "ro")

        xmin, xmax = 20000.0, 2500.0
        ymin, ymax = 1.e-4, 1.e4
        pyplot.text(xmin*.75,ymax*0.1,str(number_of_stars)+" stars\nt="+str(end_time))
        pyplot.axis([xmin, xmax, ymin, ymax])
        pyplot.savefig(name_of_the_figure)
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
        if number_of_stars > (10):
            warnings.warn("You're simulating a large number of stars with MESA. This may not be such a good idea...")
        if number_of_stars > (1000):
            raise Exception("You want to simulate with more than 1000 stars using MESA, this is not supported")
        return result
    
    def new_code(self, name_of_the_code, number_of_stars):
        if hasattr(self, name_of_the_code):
            return getattr(self, name_of_the_code)(number_of_stars)
        else:
            raise Exception("Cannot instantiate code with name '{0}'".format(name_of_the_code))

def new_code(name_of_the_code, number_of_stars):
    return InstantiateCode().new_code(name_of_the_code, number_of_stars)
        
def test_simulate_short():
    assert is_mpd_running()
    code = new_code("sse", 100)
    
    test_results_path = get_path_to_results()
    output_file = os.path.join(test_results_path, "cluster_HR_diagram.png")
    
    simulate_stellar_evolution(
            code,
            number_of_stars=100, 
            end_time = 2.0 | units.Myr,
            name_of_the_figure=output_file
    )
        
        
def new_commandline_option_parser():
    result = OptionParser(usage)
    result.add_option(
            "-c",
            "--code",
            choices=["sse", "evtwin","mesa"],
            default="sse",
            dest="code",
            metavar="CODE",
            help="CODE to use for stellar evolution"
    )
    result.add_option(
            "-n",
            "--number_of_stars",
            type="int",
            default=10,
            dest="number_of_stars",
            help="Number of stars in the cluster"
    )
    result.add_option(
            "-t",
            "--end_time",
            type="float",
            default=1000.0,
            dest="end_time",
            help="Time to evolve to in Myr"
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
            default="cluster_HR_diagram.png",
            dest="plot_filename",
            help="Name of the file to plot to"
    )
    
    return result
        
if __name__ == '__main__':
    
    parser = new_commandline_option_parser()
    (options, arguments) = parser.parse_args()
    if arguments:
        parser.error("unknown arguments '{0}'".format(arguments))
    
    code = new_code(options.code, options.number_of_stars)

    if not (options.cacheDir is None):
        print "Using cache directory: %s" % (options.cacheDir)
        code = CachedStellarEvolution(code, options.cacheDir)

    if not (options.salpeterSeed is None):
        numpy.random.seed(options.salpeterSeed)
    
    simulate_stellar_evolution(
            code,
            number_of_stars=options.number_of_stars, 
            end_time = options.end_time | units.Myr,
            name_of_the_figure=options.plot_filename
    )
