import numpy
from amuse.lab import *
from optparse import OptionParser
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel

from amuse.community.fractalcluster.interface import new_fractal_cluster_model

from prepare_figure import figure_frame, quad_frame, single_frame
from distinct_colours import get_distinct

def plot_projected_density(model, xmin=-1, xmax=1, col=0):
    pyplot.xlim(xmin, xmax)
    pyplot.ylim(xmin, xmax)
    cols = get_distinct(4)
    pyplot.scatter(model.x.value_in(nbody_system.length), model.z.value_in(nbody_system.length), c=cols[col], s=80, lw=0)

def _main(N=10, W=None): 

    f, ax1, ax2, ax3, ax4 = quad_frame("X [length]", "Y [length]")
    model = new_plummer_model(N)
    plot_projected_density(model, ax1)

    model = new_king_model(N, 9.0)
    plot_projected_density(model, ax2)

    model = new_halogen_model(N, alpha=1., beta=5., gamma=0.5)
    plot_projected_density(model, ax3)

    model = new_fractal_cluster_model(N=N, fractal_dimension=1.6)
    plot_projected_density(model, ax4)

#    pyplot.show()
    pyplot.savefig("density_distributions.eps", bbox_inches="tight")  

def plummer_model(N, x_label = 'x [length]', y_label='y [length]'):
    fig = single_frame(x_label, y_label, xsize=8, ysize=8)
    model = new_plummer_model(N)
    plot_projected_density(model, col=0)
    pyplot.savefig("plummer_model")

def king_model(N, W=9, x_label = 'x [length]', y_label='y [length]'):
    fig = single_frame(x_label, y_label, xsize=8, ysize=8)
    ax = pyplot.gca()
    model = new_king_model(N, W)
    plot_projected_density(model, col=1)
    pyplot.savefig("king_model")

def fractal_model(N, F=1.6, x_label = 'x [length]', y_label='y [length]'):
    fig = single_frame(x_label, y_label, xsize=8, ysize=8)
    ax = pyplot.gca()
    model = new_fractal_cluster_model(N=N, fractal_dimension=1.6, random_seed=42)
    plot_projected_density(model, col=2)
    pyplot.savefig("fractal_model")

def galaxy_model(N, x_label = 'x [length]', y_label='y [length]'):
    fig = single_frame(x_label, y_label, xsize=8, ysize=8)
    ax = pyplot.gca()
    model = new_halogen_model(N, alpha=1., beta=5., gamma=0.5)
    plot_projected_density(model, col=3)
    pyplot.savefig("galaxy_model")

def main(N=10): 
    numpy.random.seed(42)
    plummer_model(N)
    king_model(N)
    fractal_model(N)
    galaxy_model(N)

def new_option_parser():
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [1000]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
