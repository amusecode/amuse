import os
import numpy
from matplotlib import pyplot

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct

def read_vostok1999temperaturedata(datafile):
    t = []
    T = []
    for line in open(datafile).xreadlines():
        if '***' not in line and '#' not in line:
            sl = line.split()
            t.append(float(sl[1]))
            T.append(float(sl[3]))
    return t, T

def read_Earthorbit(datafile):
    t = []
    e = []
    a = []
    tmax = 1.e+6
    for line in open(datafile).xreadlines():
        if '#' not in line:
            sl = line.split()
            t.append(float(sl[0]))
            a.append(float(sl[1]))
            e.append(float(sl[2]))
            if t[-1]>tmax:
                break
    return t, a,  e

if __name__ in ('__main__', '__plot__'):

    try:
        amusedir = os.environ['AMUSE_DIR']
        dir = amusedir+'/examples/textbook/'
    except:
        print 'Environment variable AMUSE_DIR not set'
        dir = './'
    
    tv, T = read_vostok1999temperaturedata(dir+'vostok1999temperature.data')
    te, ae, ee = read_Earthorbit(dir+'EarthOrbit_Eps-3.data')
    
    q = []
    q0 = ae[0]*(1-ee[0]**2)
    ao = ae[0]
    T_mean =  numpy.mean(T)
    Tmin = numpy.min(T)
    Tmax = numpy.max(T)
    dT = Tmax-Tmin
    for i in range(len(te)):
        te[i] -= 0
        ae[i] = 2e+5*(ae[i]-ao)/ao 
        te[i] = te[i]/1000.
    for i in range(len(tv)):
        tv[i] = tv[i]/1000.

    from matplotlib import pyplot, rc
    
    x_label = "$a-a_0$ [AU]"
    y_label = "eccentricty"
    figure = single_frame(x_label, y_label, xsize=14, ysize=10)
    ax1 = pyplot.gca()
    cols = get_distinct(2)

    font = {'size' : 20}
    rc('font', **font)
    ax1.plot(tv, T, ls='-', c=cols[0])
    ax1.set_xlabel('time before present (kyr)')
    ax1.set_xlim(450, 0)

    # Make the y-axis label and tick labels match the line color.
    
    ax1.set_ylabel('temperature [$^o$]', color=cols[0])
    for tl in ax1.get_yticklabels():
        tl.set_color(cols[0])

    ax2 = ax1.twinx()
    ax2.plot(te, ee, cols[1])
    ax2.set_ylabel('eccentricity', color=cols[1])
    for tl in ax2.get_yticklabels():
        tl.set_color(cols[1])

    save_file = 'vostok1999temperature.png'
    pyplot.savefig(save_file)
    print "\nSaved figure in file", save_file,'\n'
    pyplot.show()

