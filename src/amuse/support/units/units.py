import numpy
from amuse.support.units.si import *
from amuse.support.units.derivedsi import *

# misc every day
minute = 60 * s
hour = 60 * minute
day = 24 * hour
yr =   named('year', 'yr', 365.242199 * day)
julianyr = named('julian yr','julianyr',365.25* day)
ms = named('meter per seconds', 'ms', m / s)
kms = named('kilometer per seconds', 'kms', km / s)

# astronomical units
angstrom = named('angstrom', 'angstrom', 1e-10*m)
AU =  named('astronomical unit', 'AU', 149597870691.0  * m)
AUd = named('AU per day','AUd', 149597870691.0  * m / day)
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

# cgs units
g = named('gram','g', 1e-3 * kg)
cm = named('centimeter','cm',0.01*m)
erg = named('energy','erg', 1e-7 * J)

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

#special unit for keys of particles
object_key = core.key_unit('object_key','key')

