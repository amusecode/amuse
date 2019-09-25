"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.lab import *

from prepare_figure import single_frame
from distinct_colours import get_distinct

def Bondi_Hoyle_Littleton_accretion_rate(Mb, vs, a0, Mdot_donor):
    a = (3.e-7 | units.MSun/units.yr)
    b = (Mb/(1|units.MSun))**2
    c = ((10|units.kms)/vs)**4
    d = ((100|units.AU)/a0)**2
    e = (Mdot_donor/(1.e-4 | units.MSun/units.yr))
    print a, b, c, d, e
    Mdot_BHL = a*b*c*d*e
    #Mdot_BHL = (3.e-7 | units.MSun/units.yr) * (Mb/(1|units.MSun))**2 ((10|units.kms)/vs)**4 * ((100|units.AU)/a0)**2 * (Mdot_donor/(1.e-4 | units.MSun/units.yr))
    return Mdot_BHL
            
    
def read_accretion_rate(filename):
    t = [] | units.yr
    m = [] | units.MSun
    n = []
    for line in open(filename).xreadlines():
        if "N accreted:" in  line:
            l = line.split()
            t.append(float(l[2])|units.Myr)
            n.append(float(l[4]))
            m.append(float(l[6])|units.MSun)
    return t, n, m
    
def plot_accretion_from_wind(filename, color):
    t, dn, dm = read_accretion_rate(filename)
    m = numpy.cumsum(dm)
    #MMoon = 3.69145063653e-08 | units.MSun
    #m /= MMoon
    m /= (1.e-9|units.MSun)

    pyplot.plot(t.value_in(units.yr), m, c=color)
#    pyplot.show()

def v_terminal_teff(temperature):
  print numpy.log10(temperature.value_in(units.K))-3.
  t4=(numpy.log10(temperature.value_in(units.K))-4.).clip(0.,1.)
  print t4
  return (30 | units.km/units.s) + ((4000 | units.km/units.s)*t4)

def main():

    """
    Mb = (2. + 1.) | units.MSun
    Mb = 2. | units.MSun
    a0 = 10. | units.AU
#    vs = 35. | units.kms
    vs = 30. | units.kms
    Mdot_donor = 0.11 |units.MSun/units.Myr
    Mdot = Bondi_Hoyle_Littleton_accretion_rate(Mb, vs, a0, Mdot_donor)
    print Mdot.in_(units.MEarth/units.yr)
    t = numpy.arange(0, 3, 0.1) | units.yr
    mdot = Mdot*t
    t += 4 | units.yr
    """
    Mb = (2. + 1.) | units.MSun
    Mb = 2. | units.MSun
    a0 = 10. | units.AU
    #vs = 30. | units.kms
#    vs = 18. | units.kms
    vs = 17.26 | units.kms
    vorb = numpy.sqrt(constants.G*Mb/a0)
    print "vorb=", vorb.in_(units.kms)
    k = vorb/vs
    m1 = 1.924785833858|units.MSun
    m2 = 1.|units.MSun
    mu = m2/(m1+m2)
    Mdot_donor  = 0.11 |units.MSun/units.Myr
    cvw = 0.
    print "k and mu:", k, mu
    Mdot = Mdot_donor  * mu**2 * k**4/(1 + k**2 + cvw**2)**(3./2.)
#    print "Mdot:", Mdot.in_(units.MEarth/units.yr)
    print "Mdot:", Mdot.in_(units.MSun/units.Myr)

    t = numpy.arange(0, 3, 0.1) | units.yr
    mdot = Mdot*t 
    t += 4.2 | units.yr

    c = get_distinct(3)

    x_label = 't [yr]'
    y_label = 'M [$10^{-9}$M$_{\odot}$]'
    figure = single_frame(x_label, y_label, logy=False, xsize=14, ysize=10)
    filename = "hydro_give_or_take.data"

    print t, mdot.value_in(units.MSun)
    mdot /= (1.e-9|units.MSun)
    pyplot.plot(t.value_in(units.yr), mdot, c=c[2])

    plot_accretion_from_wind(filename, c[0])
    filename = "hydro_give_or_take_gravity_NoG.data"
    plot_accretion_from_wind(filename, c[1])

    pyplot.savefig("hydro_accretion_from_windy_star")
#    pyplot.show()
 
if __name__ in ('__main__', '__plot__'):
    main()


