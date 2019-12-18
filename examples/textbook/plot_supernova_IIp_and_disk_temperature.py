import numpy
from matplotlib import pyplot
from amuse.lab import *

from prepare_figure import single_frame, figure_frame, set_tickmarks
from distinct_colours import get_distinct
from supernova_IIp_Lightcurve import Supernova_IIp

def read_supernova_irradiation_file(filename):
    t = []
    Tmean = []
    Tmax = []
    for line in open(filename):
        if 'Time' in line:
            sl = line.split()
            t.append(float(sl[1]))
        if 'Temperature:' in line:
            sl = line.split()
            Tmean.append(float(sl[3]))
            Tmax.append(float(sl[1]))
    return numpy.asarray(t, dtype="float"), Tmean, Tmax

def read_Earthorbit():
    t = []
    e = []
    a = []
    tmax = 1.e+6
    for line in open('Eart_Orbit_Eps-3.data'):
        if '#' not in line:
            sl = line.split()
            t.append(float(sl[0]))
            a.append(float(sl[1]))
            e.append(float(sl[2]))
            if t[-1]>tmax:
                break
    return t, a,  e

if __name__ in ('__main__', '__plot__'):
    to = 50|units.day

    t_offset = to + (((0.15|units.parsec)/(1|units.lightyear)) | units.yr)
    filename = 'SN10a.R0.15.i15.data'
    time, Tmean, Tmax = read_supernova_irradiation_file(filename)
    time += t_offset.value_in(units.day)

    t_offset = to + (((0.3|units.parsec)/(1|units.lightyear)) | units.yr)
    filename = 'SN11aof.R0.3.i45.data'
    t3pc_N7, Tmean3pc_N7, Tmax3pc_N7 = read_supernova_irradiation_file(filename)
    t3pc_N7 += t_offset.value_in(units.day)

    t_offset = to + (((0.4|units.parsec)/(1|units.lightyear)) | units.yr)
    filename = 'SN11aof.R0.4.i15.data'
    t3pc_N8, Tmean3pc_N8, Tmax3pc_N8 = read_supernova_irradiation_file(filename)
    t3pc_N8 += t_offset.value_in(units.day)

    PS1_11aof = Supernova_IIp("11aof", to)
    t = 10**numpy.arange(-2, 3., 0.01) | units.day
    L11aof = [] | units.erg/units.s
    for ti in t:
        L11aof.append(PS1_11aof.luminosity_at_time(ti))
    L11aof = numpy.log10(L11aof.value_in(units.LSun))

    PS1_10a = Supernova_IIp("10a", to)
    L10a = [] | units.erg/units.s
    for ti in t:
        L10a.append(PS1_10a.luminosity_at_time(ti))
    L10a = numpy.log10(L10a.value_in(units.LSun))
    
    from matplotlib import pyplot, rc
    x_label = "$t$ [day]"
    y_label = "L [L$_\odot$]"
    figure = single_frame(x_label, y_label, xsize=14, ysize=10)
    ax1 = pyplot.gca()
    cols = get_distinct(4)

    font = {'size' : 20}
    rc('font', **font)
    ax1.plot(t.value_in(units.day), L11aof, ls='-', c=cols[0])
    ax1.plot(t.value_in(units.day), L10a, ls='--', c=cols[0])
    ax1.set_xlabel('time [day]')
    ax1.set_ylabel('log$_{10}$(L/L$_\odot$)', color=cols[0])
    for tl in ax1.get_yticklabels():
        tl.set_color(cols[0])

    ax2 = ax1.twinx()
    ax2.plot(time, Tmean, cols[1], ls='--')
    ax2.plot(t3pc_N7, Tmean3pc_N7, cols[1])
    ax2.plot(t3pc_N8, Tmean3pc_N8, cols[1], lw=4)
    ax2.set_ylabel('mean temperature [K]', color=cols[1])
    for tl in ax2.get_yticklabels():
        tl.set_color(cols[1])

    t_cooling = [950, 1061]
    T_cooling = [1600, 800]
    ax2.plot(t_cooling, T_cooling, cols[3], lw=1)
    ax2.text(t_cooling[0]+20, T_cooling[0]-100, "cooling of 0.3 K/h", rotation=-76.5, color=cols[3])
        
    pyplot.show()
#    pyplot.savefig("supernova_IIp_and_disk_temperature")

