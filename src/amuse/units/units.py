import math




# misc every day
from amuse.units.si import *
from amuse.units.derivedsi import *
minute = named('minute', 'min', 60 * s)
hour   = named('hour',   'hr',  60 * minute)
day    = named('day',    'day', 24 * hour)
yr     = named('year',   'yr',  365.242199 * day)
julianyr = named('julian yr','julianyr',365.25* day)
ms = named('meter per seconds', 'ms', m / s)
kms = named('kilometer per seconds', 'kms', km / s)

# astronomical units
angstrom = named('angstrom', 'angstrom', 1e-10*m)
AU =  named('astronomical unit', 'AU', 149597870691.0  * m)
AUd = named('AU per day','AUd', 149597870691.0  * m / day)
parsec=named('parsec','parsec', AU / math.tan(math.pi/(180*60*60)))
kpc=named('kilo parsec','kpc',10**3 * parsec)
Mpc=named('mega parsec','Mpc',10**6 * parsec)
Gpc=named('giga parsec','Gpc',10**9 * parsec)
lightyear = named('light year', 'ly', 9460730472580.8 * km)
#lightyear = named('light year', 'ly', c*julianyr)
LSun = named('solar luminosity', 'LSun', 3.839e26 * W) 
MSun = named('solar mass', 'MSun', 1.98892e30 * kg)
MJupiter = named('jupiter mass', 'MJupiter', 1.8987e27 * kg)
RSun = named('solar radius', 'RSun', 6.955e8 * m)
RJupiter = named('jupiter radius', 'RJupiter', 71492. * km)
myr = named('million year', 'Myr', 1000000 * yr)
Myr = myr
gyr = named('giga (billion) year', 'Gyr', 1000000000 * yr)
Gyr = gyr

# cgs units
g = named('gram','g', 1e-3 * kg)
cm = named('centimeter','cm',0.01*m)
erg = named('erg','erg', 1e-7 * J)
barye = named('barye', 'Ba', 0.1*Pa)

percent = named('percent', '%', 0.01 * none)
metallicity = core.none_unit('metallicity', 'metallicity')

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
        "Massless Supernova",
        "Unknown stellar type",
        "Pre-main-sequence Star"
    ]
)

#special unit for keys of particles
object_key = core.key_unit('object_key','key')

