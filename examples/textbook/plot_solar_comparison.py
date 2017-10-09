import numpy
from amuse.lab import *
from amuse.units.optparse import OptionParser
from matplotlib import pyplot
from amuse.plot import scatter

from prepare_figure import single_frame
from distinct_colours import get_distinct

from amuse.community.seba.interface import SeBa
    
def main(t_end, mass, z, Tstar, Lstar):

    stellar_evolution_codes = [SeBa(), SSE(), MESA(), EVtwin()]
    label = ["SeBa", "SSE", "MESA", "EVtwin"]
    marker = ["o", "v", "<", ">"]

    x_label = "$(T-T_\odot)/T_\odot)$"
    y_label = "$(L-L_\odot)/L_\odot)$"
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=10)
    pyplot.xlim(-0.008, 0.004)
    color = get_distinct(6)
    pyplot.scatter([0], [0], marker="o", c=color[3], label="Sun", s=400, lw=0)

    for si in range(len(stellar_evolution_codes)):
        stellar = stellar_evolution_codes[si]
        stellar.parameters.metallicity = z

        star = Particles(1)
        star.mass = mass
        stellar.particles.add_particles(star)
        stellar.evolve_model(t_end)

        T = (stellar.particles.temperature-Tstar)/Tstar
        L = (stellar.particles.luminosity-Lstar)/Lstar
        if si==3: 
            pyplot.scatter(T, L, marker=marker[si],
                           color=color[0], label=label[si], s=300, lw=0)
        else:
            pyplot.scatter(T, L, marker=marker[si],
                           color=color[si], label=label[si], s=300, lw=0)
        stellar.stop()

    pyplot.legend(scatterpoints=1, loc=2)

    save_file = 'fig_SunComparison.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
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

