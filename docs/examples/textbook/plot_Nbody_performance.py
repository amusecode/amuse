import sys, os
import numpy
import math
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

black = '#000000'
green = '#00FF00'
red = '#FF0000'
blue = '#0000FF'
lbl = '#FF00FF'
magenta = '#00FFFF'

Integrator = ["Hermite", "MI6", "ph4", "Huayno", "ph4_GPU",
              "Gadget2", "BHTree", "Fi_3", "Bonsai"]
color = ['r', 'g', 'k', 'b', 'k', lbl, 'm', 'r', 'g', 'b', 'k', 'm']
color = get_distinct(len(color))
lstyles = ['-.', '-.', '-.', '-.', '--', '-', '-', '-', '-']
lwidth = [2, 2, 4, 2, 4, 2, 2, 2, 4]

def read_file(filename, column, keyword):
    x = []
    fptr = open(filename)
    lines = fptr.readlines()
    for line in lines:
        l = line.split()
        if l[1]==keyword:
            # print line
            x.append(float(line.split()[column]))
    fptr.close()
    return x

def main(filename=None, lim=-1):

    if filename is None: return

    try:
        amusedir = os.environ['AMUSE_DIR']
        dir = amusedir+'/examples/textbook/'
    except:
        print 'Environment variable AMUSE_DIR not set'
        dir = './'

    filename = dir+filename

    from matplotlib import pyplot, rc
    x_label = "N"
    y_label = "$t_{wall} [s]$"
    figure = single_frame(x_label, y_label, logx=True, logy=True,
                          xsize=14, ysize=10)

    npp = [12, 512]
    ntc = [1024, 2000*1024]
    tpp = map(lambda x: 0.01*x*x, npp) 
    ttc = map(lambda x: 1.e-6*(x*math.log(x)), ntc) 
    pyplot.plot(npp, tpp, c='k', ls='-', lw=4)
    pyplot.plot(ntc, ttc, c='k', ls='-', lw=4)
    pyplot.text(12, 8, '$N^2$')
    pyplot.text(2.e+3, 0.005, '$N \log (N)$')

    for ii, I in enumerate(Integrator):

        # sample line:	I= Hermite T= 1 time N= 4 M= 1.0 mass \
        #		E=  -0.250000002735 mass * length**2 * time**-2 \
        #		Q=  -0.451489790632 dE= -1.09402279371e-08 \
        #		Time= 2.04415297508 0.000992059707642

        x = read_file(filename, 6, I)
        y1 = read_file(filename, 22, I)
        y2 = read_file(filename, 23, I)
        for i in range(len(y1)):
            y1[i] *= 0.1
            y2[i] *= 0.1
        if len(x) > 0:
            pyplot.plot(x, y2, label=I, c=color[ii],
                        lw=lwidth[ii], ls=lstyles[ii])
            pyplot.scatter(x, y2, c=color[ii], lw=0, s=200)

    pyplot.legend(loc="lower right", fontsize=18)

    save_file = "Nbody_performance.png"
    pyplot.savefig(save_file)
    print "\nSaved figure in file", save_file,'\n'
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default='NbodyAMUSE.test',
                      help="output filename [NbodyAMUSE.test]")
    result.add_option("-l", dest="lim", type="float", default = -1,
                      help="boxsize")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
