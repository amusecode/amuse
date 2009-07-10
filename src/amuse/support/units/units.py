from amuse.support.units.si import *

km = k(m)

parsec = named('parsec','parsec', 3.08568025e16 * m)
AU =  named('astronomical unit', 'AU', 149597870691.0  * m)
lightyear = named('light year', 'ly', 9460730472580.8 * km)
MSun = named('solar mass', 'MSun', 1.98892e30 * kg)
minute = 60 * s
hour = 60 * minute
day = 24 * hour
yr =   named('year', 'yr', 365.2422 * day)
myr = named('million year', 'Myr', 1000000 * yr)

g = named('gram','g', 1e-3 * kg)
one = no_unit(1)

N = m**3/(kg * (s**2))
G = N(6.673e-11)

