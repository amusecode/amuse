"""
   Minimal routine for running a stellar evolution code.
"""
from amuse.lab import *
from matplotlib import pyplot 
import plot_M67Data 

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

TBSS = [6170, 6820, 6675, 7050, 6650]
LBSS = [12.259, 11.078, 12.127, 11.226, 12.892]

def single_star_evolution(M, z, model_time):

    stellar = SeBa()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=M))
    stellar.commit_particles()

    initial_luminosity = stellar.particles.luminosity
    dt = 1 | units.Myr
    time = 0 | units.Myr
    L = [] | units.LSun
    T = [] | units.K
    while stellar.particles[0].age<model_time:
        time += dt
        stellar.evolve_model(time)
        L.append(stellar.particles[0].luminosity)
        T.append(stellar.particles[0].temperature)

    final_luminosity = stellar.particles.luminosity
    print "L(t=0)=", initial_luminosity, \
          ", L (t=", stellar.particles.age, ")=", \
        final_luminosity, stellar.particles.radius
    stellar.stop()
    return L, T

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit= units.MSun,
                      dest="M", type="float",default = 1.7 | units.MSun,
                      help="stellar mass [1.0] %unit")
    result.add_option("-t", unit = units.Myr,
                      dest="model_time", type="float", 
                      default = 4000.0|units.Myr,
                      help="end time of the simulation [4.7] %unit")
    result.add_option("-z", dest="z", type="float", 
                      default = 0.02, help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    x_label = "T [$K$]"
    y_label = "L [$L_\odot$]"
    figure = single_frame(x_label, y_label, logx=False, logy=True,
                          xsize=14, ysize=10)
    color = get_distinct(6)

    L, T = single_star_evolution(M=1.4|units.MSun, z=0.04,
                                 model_time=4|units.Gyr)
    pyplot.plot(T.value_in(units.K),L.value_in(units.LSun), c=color[0])

    L, T = single_star_evolution(**o.__dict__)
    pyplot.plot(T.value_in(units.K),L.value_in(units.LSun), c=color[4])

    m67 = plot_M67Data.Cluster()
    m67.read()
    pyplot.scatter(m67.Teff, m67.L, c=color[1], lw=0, s=100)

    pyplot.scatter(TBSS, LBSS, c=color[3], lw=0, s=250)
    pyplot.scatter(TBSS[3], LBSS[3], c=color[2], lw=0, s=250)
    pyplot.xlim(9000, 4000)
    pyplot.ylim(1, 50)

    save_file = 'fig_M67withBSS_tracks.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
