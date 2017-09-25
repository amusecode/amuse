import os
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.io import read_set_from_file

def main(filename=None):
    if filename is None: return
    try:
        amusedir = os.environ['AMUSE_DIR']
    except:
        print 'Environment variable AMUSE_DIR not set'
        amusedir = '.'
    filename = amusedir+'/examples/textbook/'+filename

    pyplot.figure(figsize=(12,12))
    particles = read_set_from_file(filename, "amuse")
    for si in particles.history:
        scatter(si.x, si.y, s=100)
    xlabel("x")
    ylabel("y")

    save_file = "plot_gravity.png"
    pyplot.savefig(save_file)
    print "\nSaved figure in file", save_file,'\n'
    pyplot.show()

    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity.h5",
                      help="input filename [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


