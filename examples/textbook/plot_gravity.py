from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.io import read_set_from_file

def main(filename):
    pyplot.figure(figsize=(12,12))
    particles = read_set_from_file(filename, "amuse")
    for si in particles.history:
        scatter(si.x, si.y, s=100)
    xlabel("x")
    ylabel("y")
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


