import numpy
from amuse.lab import *
from amuse.ic.gasplummer import new_plummer_gas_model

from prepare_figure import single_frame, figure_frame
from distinct_colours import get_distinct

def binned_mean_data(r, x):
    R = numpy.arange(0, r[-1], 0.1)
    X = numpy.zeros(len(R))
    N = numpy.zeros(len(R))
    for i in range(len(R)-1):
        for j in range(len(r)):
            if r[j]>=R[i] and r[j]>=R[i+1]:
                X[i] += x[j]
                N[i] += 1.
    for i in range(len(X)):
        if X[i]>0 and N[i]>0:
            X[i] = X[i]/float(N[i])

    return R, X

def plot_ionization_fraction(pos, xion):
    r = [] | units.parsec
    x = []
    for pi, xi in zip(pos, xion):
        r.append(pi.length())
        x.append(xi)
    r, x = zip(*sorted(zip(r.value_in(units.parsec), x)))

    R, X = binned_mean_data(r, x)
    
    from matplotlib import pyplot
    x_label = "r [pc]"
    y_label = r'$\xi_{\rm ion}$'
    figure = single_frame(x_label, y_label, logx=False, logy=False,
                          xsize=14, ysize=8)
    pyplot.scatter(r, x, c=get_distinct(1), lw=0, s=100)
    pyplot.plot(R, X, c=get_distinct(2)[1], lw=2)
    pyplot.xlim(0, 6)
    pyplot.ylim(-0.04, 1.19)
    pyplot.savefig("fig_ionization_of_GMC")
    pyplot.show()

###BOOKLISTSTART1###
def generate_ism_initial_conditions(N, boxsize):
    converter = nbody_system.nbody_to_si(10|units.MSun, 3|units.parsec)
    ism = new_plummer_gas_model(N, converter)
    ism.flux = 0. | units.s**-1
    ism.xion = 0.0

    hydro = Fi(converter)
    hydro.gas_particles.add_particles(ism)
    hydro.evolve_model(1|units.hour)
    hydro.gas_particles.new_channel_to(ism).copy()
    hydro.stop()
    ism = ism.select(lambda r: r.length() < 0.5*boxsize,["position"])
    print "Max density:", ism.rho.max().in_(units.MSun/units.parsec**3), \
          ism.rho.max().in_(units.amu/units.cm**3)
    return ism
###BOOKLISTSTOP1###

###BOOKLISTSTART###
def main(N, Lstar, boxsize, t_end):
    ism = generate_ism_initial_conditions(N, boxsize)

    source = Particle()
    source.position = (0, 0, 0) |units.parsec
    source.flux = Lstar/(20. | units.eV)
    source.rho = ism.rho.max()
    source.xion = ism.xion.max()
    source.u = (9.|units.kms)**2

    radiative = SimpleX()
    radiative.parameters.box_size = boxsize    
    radiative.particles.add_particle(source)
    radiative.particles.add_particles(ism)

    radiative.evolve_model(t_end)
    print "min ionization:", radiative.particles.xion.min()
    print "average ionization:", radiative.particles.xion.mean()
    print "max ionization:", radiative.particles.xion.max()
    plot_ionization_fraction(radiative.particles.position,
                             radiative.particles.xion)
    radiative.stop()
###BOOKLISTSTOP###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default = 10000,
                      help="number of gas particles [%default]")
    result.add_option("-t", unit=units.Myr,
                      type="float",
                      dest="t_end", default = 1.0|units.Myr,
                      help="radiation time [%default]")
    result.add_option("-L", unit=units.LSun,
                      type="float",
                      dest="Lstar", default = 100|units.LSun,
                      help="luminosity of ionizing source [%default]")
    result.add_option("-d", unit=units.parsec,
                      dest="boxsize", default = 10|units.parsec,
                      help="size of the density box [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    numpy.random.seed(12345)
    main(**o.__dict__)

