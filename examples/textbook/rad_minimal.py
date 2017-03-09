from amuse.lab import *
from amuse.ic.gasplummer import new_plummer_gas_model

from prepare_figure import single_frame, figure_frame
from distinct_colours import get_distinct

def plot_ionization_fraction(pos, xion):
    r = [] | units.parsec
    x = []
    for pi, xi in zip(pos, xion):
        r.append(pi.length())
        x.append(xi)
    r, x = zip(*sorted(zip(r.value_in(units.parsec), x)))
    
    from matplotlib import pyplot
    x_label = "r [pc]"
    y_label = r'$\xi_{\rm ion}$'
    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=14, ysize=8)
    pyplot.scatter(r, x, c=get_distinct(1), lw=0, s=100)
    pyplot.xlim(0, 10)
    pyplot.ylim(-0.04, 1.19)
    #pyplot.savefig("fig_ionization_of_GMC")
    pyplot.show()

def main(N, Lstar, boxsize, rho, t_end):

    source=Particle()
    source.position = (0, 0, 0) |units.parsec
    source.flux = Lstar/(20. | units.eV)
    source.rho = rho
    source.xion = 0.0
    source.u = (9. |units.kms)**2

    converter=nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
    ism = new_plummer_gas_model(N, converter)
    ism.rho = rho
    ism.u = source.u
    ism.flux = 0. | units.s**-1
    ism.xion = source.xion
    ism = ism.select(lambda r: r.length()<0.5*boxsize,["position"])

    radiative = SimpleX()
    radiative.parameters.box_size=1.001*boxsize    
    radiative.parameters.timestep=0.001 | units.Myr

    radiative.particles.add_particle(source)
    radiative.particles.add_particles(ism)

    radiative.evolve_model(t_end)
    print "min ionization:", radiative.particles.xion.min()
    print "average ionization:", radiative.particles.xion.mean()
    print "max ionization:", radiative.particles.xion.max()

    plot_ionization_fraction(radiative.particles.position, radiative.particles.xion)
    radiative.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default = 10000,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", default = 0.1|units.Myr,
                      help="radiation time [%default]")
    result.add_option("-L", unit=units.LSun,
                      dest="Lstar", default = 100|units.LSun,
                      help="luminosity of ionizing source [%default]")
    result.add_option("-p", unit=units.amu/units.cm**3, 
                      dest="rho", default = 1|units.amu/units.cm**3,
                      help="interstellar density [%default] amu/cm^3")
    result.add_option("-d", unit=units.parsec,
                      dest="boxsize", default = 100|units.parsec,
                      help="size of the density box [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

