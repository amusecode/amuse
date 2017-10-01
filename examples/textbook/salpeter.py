from amuse.lab import new_powerlaw_mass_distribution

def generate_power_law_mass_function(N, Mmin, Mmax, ximf):
    masses = new_powerlaw_mass_distribution(N, Mmin, Mmax, ximf)
    plot_mass_function(masses, ximf)

###BOOKLISTSTART###
import numpy
import math
from amuse.units import units
from matplotlib import pyplot
from prepare_figure import figure_frame, get_distinct

def plot_mass_function(masses, ximf):
    Mmin = masses.min()
    Mmax = masses.max()
    lm = math.log10(0.5*Mmin.value_in(units.MSun))
    lM = math.log10(1.5*Mmax.value_in(units.MSun))
    bins = 10**numpy.linspace(lm, lM, 51)
    Nbin, bin_edges= numpy.histogram(masses.value_in(units.MSun), bins=bins)
    y = Nbin / (bin_edges[1:] - bin_edges[:-1])
    x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    for i in range(len(y)):
        y[i] = max(y[i], 1.e-10)

    fig, ax = figure_frame("M$_\odot$", "N", xsize=12, ysize=8)
    colors = get_distinct(2)
    pyplot.scatter(x, y, s=100, c=colors[0], lw=0)

    c = ((Mmax.value_in(units.MSun)**(ximf+1)) \
         	- (Mmin.value_in(units.MSun)**(ximf+1))) / (ximf+1)
    pyplot.plot(x, len(masses)/c * (x**ximf), c=colors[1])
    pyplot.loglog()
    
    save_file = "salpeter.png"
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
###BOOKLISTSTOP###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [10]")
    result.add_option("-m", unit=units.MSun, 
                      dest="Mmin", type="float",default = 1|units.MSun,
                      help="minimum mass of the mass function [0.1] %unit")
    result.add_option("-M", unit=units.MSun, 
                      dest="Mmax", type="float",default = 100|units.MSun,
                      help="maximum mass of the mass function [100] %unit")
    result.add_option("-x", 
                      dest="ximf", type="float",default = -2.35,
                      help="mass function slope [-2.35]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    numpy.random.seed(31415)
    generate_power_law_mass_function(**o.__dict__)

