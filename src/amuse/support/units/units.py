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

N = named('Newton', 'N', kg * m /s**2)
J = named('joule','J', kg * m **2  * s ** -2)
W = named('watt', 'W', J / s)
F = named('farad','F', s**4*A**2*m**-2*kg**-1)
C = named('coulomb','C', A*s)
V = named('volt','V', J/C)

amu=named('atomic mass unit', 'amu',1.66053886*10**-27 * kg)
e=named('electron charge','e',1.6021765314e-19 * C)
eV=named('electron volt','eV', e*V)

# cgs
g = named('gram','g', 1e-3 * kg)
cm=named('centimeter','cm',0.01*m)
erg = named('energy','erg', 1e-7 * J)

# physical constants
G = 6.673e-11 | m**3/kg/s**2
kboltz = 1.3806503 * 10**-23 | m**2 * kg / s**2 / K
c=299792458. | m / s
h=6.6260689633e-34 | J * s
hbar=h/2/numpy.pi
eps0=8.854187817e-12 | F/m


# astronomical units
AU =  named('astronomical unit', 'AU', 149597870691.0  * m)
AUd = named('AU per day','AUd', 149597870691.0  * m / day)
#parsec = named('parsec','parsec', 3.08568025e16 * m)
parsec=named('parsec','parsec', 3.08568025e16 * m)
kpc=named('kilo parsec','kpc',10**3 * parsec)
Mpc=named('mega parsec','Mpc',10**6 * parsec)
lightyear = named('light year', 'ly', 9460730472580.8 * km)
#lightyear = named('light year', 'ly', c*julianyr)
LSun = named('solar luminosity', 'LSun', 3.839e26 * W) 
MSun = named('solar mass', 'MSun', 1.98892e30 * kg)
RSun = named('solar radius', 'RSun', 6.955e8 * m)
myr = named('million year', 'Myr', 1000000 * yr)
Myr = myr

Pa = named('Pascal', 'Pa', N / (m ** 2))
weber = named('weber', 'Wb', kg * m ** 2 * s ** -2 * A ** -1) 
tesla = named('tesla', 'T', weber / (m ** 2))

percentage = core.none_unit('percentage', '%')
string = core.string_unit('string', 'string')

stellar_type = core.enumeration_unit(
    'stellar_type',
    'stellar_type',
    None,
    [
        "deeply or fully convective low mass MS star",
        "Main Sequence star",
        "Hertzsprung Gap",
        "First Giant Branch",
        "Core Helium Burning",
        "First Asymptotic Giant Branch",
        "Second Asymptotic Giant Branch",
        "Main Sequence Naked Helium star",
        "Hertzsprung Gap Naked Helium star",
        "Giant Branch Naked Helium star",
        "Helium White Dwarf",
        "Carbon/Oxygen White Dwarf",
        "Oxygen/Neon White Dwarf",
        "Neutron Star",
        "Black Hole",
        "Massless Supernova"
    ]
)
