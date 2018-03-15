#!/usr/bin/python

from numpy import *
from amuse.lab import *

#a XiTau = 263.973|units.RSun = 1.2272459225 AU
#HD a= 171.109773416 RSun = quantity<0.795511639712 AU>

def orbital_separation_after_mass_transfer(a, md, mdf, ma, maf, nj) :
    mt = ma+md
    mtf = maf+mdf
    fa = ((mdf*maf)/(md*ma))**(-2.) * (mtf/mt)**(2.*nj+1.)
#    print "double_star.orbital_separation_after_mass_transfer()", fa
    return a*fa


def alpha_lambda_due_to_mass_ejection(ai, af, dM, m1, dm1, m2, dm2):
    al = -1 * ((m1+m1)*dM/ai) / (2. * ( (m1+dm1)*(m2+dm2)/ai - m1*m2/af))
    return al

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="mdonor",type="float",default=5.50|units.MSun)
    result.add_option("-m", unit=units.MSun,
                      dest="macc",type="float",default=6.30|units.MSun)
    result.add_option("--dMd", unit=units.MSun,
                      dest="dmd",type="float",default=-0.06|units.MSun)
    result.add_option("--dma", unit=units.MSun,
                      dest="dma",type="float",default=0.02|units.MSun)
    result.add_option("-a", unit=units.AU,
                      dest="a",type="float",default=1.2272459225|units.AU)
    result.add_option("--dt", unit=units.day,
                      dest="dt",type="float",default=1400.|units.day)
    result.add_option("--nj", 
                      dest="nj",type="float",default=3)
    return result

def main(dt, a, mdonor, dmd, macc, dma, nj):
    mdonor_new = mdonor+dmd
    macc_new = macc+dma
    a_new = orbital_separation_after_mass_transfer(a, mdonor, mdonor_new, 
                                                   macc, macc_new, nj)
    print "a_old=", a.as_quantity_in(units.AU)
    print "a_new=", a_new.as_quantity_in(units.AU)

if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)

