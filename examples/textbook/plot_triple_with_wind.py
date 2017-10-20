import numpy
from amuse.lab import *

from prepare_figure import single_frame
from distinct_colours import get_distinct
from matplotlib import pyplot

def read_triple_data(filename):
    t = [] 
    ain = [] 
    aout = []
    ein = []
    eout = []
    a0in = 0
    a0out = 0
    for line in open(filename).xreadlines():
        if "Triple" in  line:
            l = line.split()
            ti = float(l[3])
            if ti<=0:
                a0in = float(l[10])
                a0out = float(l[16])
                e0in = float(l[12])
                e0out = float(l[18])
            if ti>=4:
                t.append(float(l[3]))
                ain.append(float(l[10])/a0in)
                ein.append(float(l[12])/e0in)
                aout.append(float(l[16])/a0out)
                eout.append(float(l[18])/e0out)
    return t, ain, ein, aout, eout

filename = "evolve_triple_with_wind.data"
#filename = "a"
t, ain, ein, aout, eout = read_triple_data(filename)
print len(ain), len(aout)

x_label = "$a/a_{0}$"
y_label = "$e/e_{0}$"
fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=14, ysize=12)
color = get_distinct(2)

pyplot.plot(ain, ein, c=color[0], label= 'inner')
pyplot.plot(aout, eout, c=color[1], label= 'outer')
pyplot.legend(loc="upper left", ncol=1, shadow=False, fontsize=24)


#pyplot.show()
pyplot.savefig("evolve_triple_with_wind")
