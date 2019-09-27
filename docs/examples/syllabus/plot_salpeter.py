"""
   Example AMUSE sciprt for plotting the Salpeter mass function
"""
import numpy 
from matplotlib import pyplot
from amuse.lab import new_salpeter_mass_distribution, units

def main(N, m, M, ximf):
    masses = new_salpeter_mass_distribution(N, m, M, ximf)
    lm = numpy.log10(m.value_in(units.MSun))
    lM = numpy.log10(M.value_in(units.MSun))
    bins = 10**numpy.linspace(lm, lM, 50)
    Nbin, bin_edges= numpy.histogram(masses.value_in(units.MSun), bins=bins)
    bin_sizes = bin_edges[1:] - bin_edges[:-1]
    y = Nbin / bin_sizes
    x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    for i in range(len(y)):
        y[i] = max(y[i], 1.e-10)
    pyplot.scatter(x, y)
    
    c = ((M.value_in(units.MSun)**(ximf+1)) - (m.value_in(units.MSun)**(ximf+1))) / (ximf+1)
    pyplot.plot(x, N/ c * (x**ximf))
    pyplot.loglog()
    pyplot.xlabel('$M [M_\odot]$')
    pyplot.ylabel('N')
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
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

