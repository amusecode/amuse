""" 
   Example AMUSE script to generate a Plummer sphere and plot the results.
"""

###BOOKLISTSTART###
from matplotlib.pyplot import show, xlim, ylim, figure
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import new_plummer_model

def main(N=10): 
    figure(figsize=(5,5))
    bodies = new_plummer_model(N)
    scatter(bodies.x, bodies.y)
    xlim(-1, 1)
    ylim(-1, 1)
    xlabel("X")
    ylabel("Y")
    show()
###BOOKLISTSTOP###
    
def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default=1000,
                      help="number of stars [1000]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
