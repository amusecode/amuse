## The original code is in: /home/spz/Latex/papers/2012/PZMcMvE/data

import sys
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

#Integrator = ["Hermite", "MI6", "Huayno", "ph4", "ph4_GPU", "PhiGPU", "BHTree", "Fi_3", "Gadget2", "Bonsai"]
Integrator = ["Hermite", "MI6", "ph4", "Huayno", "ph4_GPU", "Gadget2", "BHTree", "Fi_3", "Bonsai"]
color = ['r', 'g', 'k', 'b', 'k', lbl, 'm', 'r', 'g', 'b', 'k', 'm']
color = get_distinct(len(color))
lstyles = ['-.', '-.', '-.', '-.', '--', '-', '-', '-', '-']
lwidth = [2, 2, 4, 2, 4, 2, 2, 2, 4]
#color='#112233'
def read_file(filename, column, keyword):
    x = []
    fptr = open(filename)
    lines = fptr.readlines()
    for line in lines:
        l = line.split()
        if l[1]==keyword:
#        if line.find(keyword)>=0:
            print line
            x.append(float(line.split()[column]))
    fptr.close()
    return x

def main(filename = "NbodyAMUSE.test", lim=-1):

#I= Hermite T= 1 time N= 4 M= 1.0 mass E=  -0.250000002735 mass * length**2 * time**-2 Q=  -0.451489790632 dE= -1.09402279371e-08 Time= 2.04415297508 0.000992059707642

    from matplotlib import pyplot, rc
    x_label = "N"
    y_label = "$t_{wall} [s]$"
    figure = single_frame(x_label, y_label, logx=True, logy=True, xsize=12, ysize=10)

    npp = [12, 512]
    ntc = [1024, 2000*1024]
    tpp = map(lambda x: 0.01*x*x, npp) 
    ttc = map(lambda x: 1.e-6*(x*math.log(x)), ntc) 
    pyplot.plot(npp, tpp, c='k', ls='-', lw=4)
    pyplot.plot(ntc, ttc, c='k', ls='-', lw=4)
    pyplot.text(12, 8, '$N^2$')
    pyplot.text(2.e+3, 0.005, '$N \log (N)$')

    for ii, I in enumerate(Integrator):
        x = read_file(filename, 6, I)
        y1 = read_file(filename, 22, I)
        y2 = read_file(filename, 23, I)
        for i in range(len(y1)):
            y1[i] *= 0.1
            y2[i] *= 0.1
#        pyplot.plot(x, y1, label=I, c=color[ii], lw=1,ls=lstyles[ii])
#        pyplot.scatter(x, y1, c=color[ii])
        if len(x)>0:
            pyplot.plot(x, y2, label=I, c=color[ii], lw=lwidth[ii], ls=lstyles[ii])
            pyplot.scatter(x, y2, c=color[ii], lw=0, s=200)
#        pyplot.plot(x, y1, label=I, c=color[ii], lw=2,ls=lstyles[ii])
#        pyplot.plot(x, y2, c=color[ii], lw=3.0,ls=lstyles[ii])
#    pyplot.legend(loc="upper left")
######    pyplot.legend(loc="lower right")
#    pyplot.xlabel("N")
#    pyplot.ylabel("t$_{wall}$ [s]")
#    pyplot.loglog()
#    pyplot.show()
    pyplot.savefig("Nbody_performance")

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "NbodyAMUSE.test",
                      help="output filename [NbodyAMUSE.test]")
    result.add_option("-l", dest="lim", type="float", default = -1,
                      help="boxsize")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
