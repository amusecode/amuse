"""
Units supported in AMUSE
"""

import numpy
from . import constants
from . import physical_constants
from . import astronomical_constants
from . import quantities

# The two imports below are to explicitly expose everything directly used in
# this module.
from .si import (named, s, m, kg, none, core, no_unit,)
from .derivedsi import (km, V, W, J, Pa, rad,)

# importing * from si and derivedsi since we do want all of the units in these
# modules, and we don't want to repeat everything here.
from .si import *
from .derivedsi import *
from .stellar_types import stellar_type

# misc every day
minute = named("minute", "min", 60.0 * s)
hour = named("hour", "hr", 60.0 * minute)
day = named("day", "day", 24.0 * hour)
yr = named("year", "yr", 365.242199 * day)
julianyr = named("julian yr", "julianyr", 365.25 * day)
ms = named("meter per seconds", "ms", m / s)
kms = named("kilometer per seconds", "kms", km / s)

# units based on measured quantities
e = named("electron charge", "e", constants.elementary_charge.as_unit())
eV = named("electron volt", "eV", e * V)
MeV = named("mega electron volt", "MeV", 1e6 * eV)
GeV = named("giga electron volt", "GeV", 1e9 * eV)
E_h = named("hartree energy", "E_h", constants.Hartree_energy.as_unit())
amu = named("atomic mass unit", "amu", constants.u.as_unit())
Ry = named(
    "rydberg unit",
    "Ry",
    (constants.Rydberg_constant * constants.h * constants.c)
    .as_quantity_in(eV)
    .as_unit(),
)

# astronomical units
angstrom = named("angstrom", "angstrom", 1e-10 * m)
au = named("astronomical unit", "au", astronomical_constants.au.as_unit())
aud = named("au per day", "aud", astronomical_constants.au.as_unit() / day)
parsec = named("parsec", "parsec", astronomical_constants.parsec.as_unit())
kpc = named("kilo parsec", "kpc", 10**3 * parsec)
Mpc = named("mega parsec", "Mpc", 10**6 * parsec)
Gpc = named("giga parsec", "Gpc", 10**9 * parsec)
lightyear = named("light year", "ly", (physical_constants.c * julianyr).as_unit())
LSun = named("solar luminosity", "LSun", astronomical_constants.Lsun.as_unit())
MSun = named("solar mass", "MSun", astronomical_constants.Msun.as_unit())
RSun = named("solar radius", "RSun", astronomical_constants.Rsun.as_unit())
MJupiter = named("jupiter mass", "MJupiter", astronomical_constants.Mjupiter.as_unit())
RJupiter = named("jupiter radius", "RJupiter", astronomical_constants.Rjupiter.as_unit())
MEarth = named("earth mass", "MEarth", astronomical_constants.Mearth.as_unit())
REarth = named("earth radius", "REarth", astronomical_constants.Rearth.as_unit())
kyr = named("kilo year", "kyr", 1000 * yr)
myr = named("million year", "Myr", 1000000 * yr)
gyr = named("giga (billion) year", "Gyr", 1000000000 * yr)

# alternatives, for compatibility
AU = au
AUd = aud
Myr = myr
Gyr = gyr
pc = parsec
Lsun = LSun
Msun = MSun
Rsun = RSun
Mjupiter = MJupiter
Rjupiter = RJupiter
Mearth = MEarth
Rearth = REarth

# cgs units
g = named("gram", "g", 1e-3 * kg)
cm = named("centimeter", "cm", 0.01 * m)
erg = named("erg", "erg", 1e-7 * J)
barye = named("barye", "Ba", 0.1 * Pa)

# imperial distance units
inch = named("inch", "in", 0.0254 * m)
foot = named("foot", "ft", 0.3048 * m)
mile = named("mile", "mi", 1609.344 * m)

percent = named("percent", "%", 0.01 * none)
metallicity = core.none_unit("metallicity", "metallicity")

string = core.string_unit("string", "string")

# special unit for keys of particles
object_key = core.key_unit("object_key", "key")

# angles
pi = numpy.pi | rad
rev = named("revolutions", "rev", (2 * numpy.pi) * rad)
deg = named("degree", "deg", (numpy.pi / 180) * rad)
arcmin = named("arcminutes", "arcmin", (1.0 / 60) * deg)
arcsec = named("arcseconds", "arcsec", (1.0 / 3600) * deg)
