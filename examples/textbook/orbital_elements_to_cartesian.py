from amuse.lab import *
from math import pi, sin, cos, sqrt, atan2
import numpy

#Solve Kepler equation by iterating: M = E - e sin E
#Lit.: Sterne, T.E., 1960, An introduction to Celestial Mechanics, p. 13-14
def eccentric_anomaly(mean_anomaly, e) :
    ecc_anomaly = mean_anomaly + 2*pi*(e * sin(mean_anomaly) + 0.5*e*e*sin(2*mean_anomaly))
    m = ecc_anomaly - e*sin(ecc_anomaly)
    de = (mean_anomaly-m) / (1 - e*cos(ecc_anomaly))
    ecc_anomaly += de;
    while de >= 0.001 :
        m = ecc_anomaly - e*sin(ecc_anomaly)
        de = (mean_anomaly-m) / (1 - e*cos(ecc_anomaly))
        ecc_anomaly += de
    return ecc_anomaly

def orbital_elements_to_pos_and_vel(time, a, ecc, inc, omra, omega, tp, P, Mbh, mstar):

    mu = constants.G*(Mbh+mstar)   
    MA = 2.*pi*(time-tp)/P
    while (MA<0) and (MA>2.*pi): 
        if MA<0: time = time+P
        else: time = time-P
        MA = 2.*pi*(time-tp)/P
  
    EA = eccentric_anomaly(MA,ecc)  # eccentric anomaly
    # true anomaly in the correct quadrant     
    ta    = 2.*atan2(sqrt(1.+ecc)*sin(EA/2.),sqrt(1.-ecc)*cos(EA/2.))
    radius = a * (1. - ecc*cos(EA)) # radius from EA and ecc

    r = [] | units.AU # Cartesian position
    r.append(radius*(cos(omra)*cos(omega+ta) - sin(omra)*sin(omega+ta)*cos(inc)))
    r.append(radius*(sin(omra)*cos(omega+ta) + cos(omra)*sin(omega+ta)*cos(inc)))
    r.append(radius*(sin(inc)*sin(omega+ta)))

    h = (mu*a*(1. - ecc*ecc)).sqrt()
    pp = a*(1-ecc*ecc)

    v = [] | units.kms # Cartesian velocity
    v.append(r.x*h*ecc/radius/pp*sin(ta) - h/radius * ( cos(omra)*sin(omega+ta) +sin(omra)*cos(omega+ta)*cos(inc)))
    v.append(r.y*h*ecc/radius/pp*sin(ta) - h/radius * ( sin(omra)*sin(omega+ta) -cos(omra)*cos(omega+ta)*cos(inc)))
    v.append(r.z*h*ecc/radius/pp*sin(ta) + h/radius*sin(inc)*cos(omega+ta))

    return r, v

def main(T, a, e, i, o, O, t, P, M, m):
    T = T |units.yr
    a = a |units.AU
    t = t | units.yr
    P = P | units.yr
    M = M |units.MSun
    m = m |units.MSun
    i *= pi/180.
    o *= pi/180.
    O *= pi/180.
    r, v = orbital_elements_to_pos_and_vel(T, a, e, i, o, O, t, P, M, m)
    print "r=", r.in_(units.AU), "v=", v.in_(units.kms)

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    # data for S2 from 2009ApJ...692.1075G
    result.add_option("-T", dest="T",type="float",default=0)
    result.add_option("-a", dest="a",type="float",default=1042.5)
    result.add_option("-e", dest="e",type="float",default=0.88)
    result.add_option("-i", dest="i",type="float",default=135.25)
    result.add_option("-o", dest="o",type="float",default=225.39)
    result.add_option("-O", dest="O",type="float",default=63.56)
    result.add_option("-t", dest="t",type="float",default=2002.32)
    result.add_option("-P", dest="P",type="float",default=15.8)
    result.add_option("-M", dest="M",type="float",default=4.45e+6)
    result.add_option("-m", dest="m",type="float",default=19.5)
    return result

if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    main(**options.__dict__)

