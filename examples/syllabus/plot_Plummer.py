from amuse.lab import *
from optparse import OptionParser
from matplotlib.pyplot import show, xlim, ylim, figure
from amuse.plot import scatter, xlabel, ylabel

from amuse.community.fractalcluster.interface import new_fractal_cluster_model

def main(N=10): 
    figure(figsize=(10,10))
    bodies = new_plummer_model(N)
    scatter(bodies.x, bodies.y)
    xlim(-1, 1)
    ylim(-1, 1)
    xlabel("X")
    ylabel("Y")
    show()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [1000]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
