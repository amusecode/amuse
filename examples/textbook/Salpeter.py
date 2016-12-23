import numpy 
from matplotlib import pyplot
from amuse.lab import *

def generate_and_bin_Salpeter(N, m, M, ximf):
    numpy.random.seed(31415)
    masses = new_salpeter_mass_distribution(N, m, M, ximf)
    m_lower = numpy.log10(0.5*m.value_in(units.MSun))
    M_upper = numpy.log10(1.5*M.value_in(units.MSun))
    bins = 10**numpy.linspace(m_lower, M_upper, 50)
    Nbin, bin_edges= numpy.histogram(masses.value_in(units.MSun), bins=bins)
    bin_sizes = bin_edges[1:] - bin_edges[:-1]
    y_binned = Nbin / bin_sizes
    M_binned = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    for i in range(len(y_binned)):
        y_binned[i] = max(y_binned[i], 1.e-10)

    return M, M_binned, y_binned

def plot_Salpeter(N, m, M, ximf):
    M, x, y = generate_and_bin_Salpeter(N, m, M, ximf)
    
    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(16, 10))
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on() 
    ax.locator_params(nbins=3)
    pyplot.xlabel("$M [M_\odot]$")
    pyplot.ylabel("N")

    pyplot.scatter(x, y, s=80, c='b')
    c = ((M.value_in(units.MSun)**(ximf+1)) - (m.value_in(units.MSun)**(ximf+1))) / (ximf+1)
    pyplot.plot(x, N/ c * (x**ximf), c='r')
    pyplot.loglog()
    pyplot.show()

from prepare_figure import figure_frame, single_frame
from distinct_colours import get_distinct

def plot_Salpeter_nice(N, m, M, ximf, output_filename):
    M, x, y = generate_and_bin_Salpeter(N, m, M, ximf)
    
    x_label = "$M [M_\odot]$"
    y_label = "N"
    fig, ax = figure_frame(x_label, y_label, xsize=12, ysize=8)
    cols = get_distinct(2)

    pyplot.scatter(x, y, s=80, c=cols[0])
    c = ((M.value_in(units.MSun)**(ximf+1)) - (m.value_in(units.MSun)**(ximf+1))) / (ximf+1)
    pyplot.plot(x, N/ c * (x**ximf), c=cols[1])
    pyplot.loglog()

    if output_filename:
        pyplot.savefig(output_filename)
    else:
        pyplot.show()
        
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [10]")
    result.add_option("-m", unit=units.MSun, 
                      dest="m", type="float",default = 1|units.MSun,
                      help="minimum mass of the mass function [0.1] %unit")
    result.add_option("-M", unit=units.MSun, 
                      dest="M", type="float",default = 100|units.MSun,
                      help="maximum mass of the mass function [100] %unit")
    result.add_option("-x", 
                      dest="ximf", type="float",default = -2.35,
                      help="mass function slope [-2.35]")
    result.add_option("-o", 
                      dest="output_filename", default ="salpeter",
                      help="output filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    plot_Salpeter_nice(**o.__dict__)

