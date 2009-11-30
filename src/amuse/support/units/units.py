from amuse.support.units.si import *

km = k(m)

parsec = named('parsec','parsec', 3.08568025e16 * m)
AU =  named('astronomical unit', 'AU', 149597870691.0  * m)
lightyear = named('light year', 'ly', 9460730472580.8 * km)
MSun = named('solar mass', 'MSun', 1.98892e30 * kg)
RSun = named('solar radius', 'RSun', 6.955e8 * m)
minute = 60 * s
hour = 60 * minute
day = 24 * hour
yr =   named('year', 'yr', 365.2422 * day)
myr = named('million year', 'Myr', 1000000 * yr)
Myr = myr
AUd = named('AU per day','AUd', 149597870691.0  * m / 365 / day)
ms = named('meter per seconds', 'ms', m / s)
g = named('gram','g', 1e-3 * kg)
one = 1 | none


G = 6.673e-11 | m**3/kg/s**2
N = named('Newton', 'N', kg / m /s**2)
J = named('joule','J', kg * m **2  * s ** -2)
erg = named('energy','erg', 1e-7 * J)
W = named('watt', 'W', J / s)

LSun = named('solar luminotisity', 'LSun', 3.839e26 * W) 

Pa = named('Pascal', 'Pa', N / (m ** 2))
weber = named('weber', 'Wb', kg * m ** 2 * s ** -2 * A ** -1) 
tesla = named('tesla', 'T', weber / (m ** 2))

percentage = core.none_unit('percentage', '%')
string = core.string_unit('string', 'string')

stellar_type = core.enumeration_unit(
    'stellar type',
    'stellar type',
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
