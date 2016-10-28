import numpy
from amuse.lab import *
from amuse.units.optparse import OptionParser
from matplotlib import pyplot
from amuse.plot import scatter

from amuse.community.seba.interface import SeBa
    
def main(t_end, mass, z, Tstar, Lstar):

    stellar_evolution_codes = [SeBa(), SSE(), MESA(), EVtwin()]
#    stellar_evolution_codes = [SeBa(), SSE()]
    label = ["SeBa", "SSE", "MESA", "EVtwin"]
    marker = ["o", "v", "<", ">"]
    color = ["r","b", "k", "m"]

    from matplotlib import pyplot, rc
    figure = pyplot.figure(figsize=(10,10))
    font = {'size' : 20}
    rc('font', **font)
    plot = figure.add_subplot(1,1,1)
    pyplot.scatter([0], [0], marker="o", c="y", label="Sun", s=200)

    for si in range(len(stellar_evolution_codes)):
        stellar = stellar_evolution_codes[si]
        stellar.parameters.metallicity = z

        star = Particles(1)
        star.mass = mass
        stellar.particles.add_particles(star)
        stellar.evolve_model(t_end)

        T = (stellar.particles.temperature-Tstar)/Tstar
        L = (stellar.particles.luminosity-Lstar)/Lstar
        pyplot.scatter(T, L, marker=marker[si], color=color[si], label=label[si], s=200)
        stellar.stop()

    pyplot.legend(scatterpoints=1, loc=2)

    pyplot.xlabel("$(T-T_\odot)/T_\odot)$")
    pyplot.ylabel("$(L-L_\odot)/L_\odot)$")
    pyplot.show()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-T", unit=units.K,
                      dest="Tstar", type="float",default = 5778 |units.K,
                      help="stellar temperature [%defailt]")
    result.add_option("-L", unit=units.LSun,
                      dest="Lstar", type="float",default = 1 |units.LSun,
                      help="stellar luminosity [%defailt]")
    result.add_option("-m", unit=units.MSun,
                      dest="mass", type="float",default = 1.0 |units.MSun,
                      help="stellar mass [%defailt]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 4.6 |units.Gyr,
                      help="end time of the simulation [%defailt]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [%defailt]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

