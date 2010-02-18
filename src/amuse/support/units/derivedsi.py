import numpy
from amuse.support.units.si import *

# handy definitions
one = 1 | none
km = k(m)

# misc every day
minute = 60 * s
hour = 60 * minute
day = 24 * hour
yr =   named('year', 'yr', 365.242199 * day)
julianyr = named('julian yr','julianyr',365.25* day)
ms = named('meter per seconds', 'ms', m / s)
kms = named('kilometer per seconds', 'kms', km / s)
Hz = named('Hertz', 'Hz', 1/s)
MHz = named('MegaHertz', 'MHz', 1e6*Hz)

N = named('Newton', 'N', kg * m /s**2)
J = named('joule','J', kg * m **2  * s ** -2)
W = named('watt', 'W', J / s)
F = named('farad','F', s**4*A**2*m**(-2)*kg**(-1))
C = named('coulomb','C', A*s)
V = named('volt','V', J/C)
T = named('Tesla','T', kg/A/s/s)
ohm = named('Ohm','ohm', V/A)
Wb = named('Weber','Wb', V*s)
sr = named('Steradian','sr',m**2/m**2)

#physical constants...
#amu=named('atomic mass unit', 'amu',1.66053886*10**-27 * kg)
#u = amu
e=named('electron charge','e',1.6021765314e-19 * C)
eV=named('electron volt','eV', e*V)
MeV=named('mega electron volt','eV', 1e6*eV)
GeV=named('giga electron volt','GeV', 1e9*eV)
E_h = named('hartree_energy','E_h', 4.359744e-18 * J)
S = named('Siemens', 'S', A/V)

# cgs
g = named('gram','g', 1e-3 * kg)
cm = named('centimeter','cm',0.01*m)
erg = named('energy','erg', 1e-7 * J)
