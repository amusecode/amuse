"""
Series of astronomical constants, which are in turn used by the units of the
same names.
"""

import numpy as np
from amuse.units.si import m, kg
from amuse.units.derivedsi import W, km


au = 149597870691.0 | m
parsec = au / np.tan(np.pi / (180 * 60 * 60))
Lsun = 3.839e26 | W
Msun = 1.98892e30 | kg
Rsun = 6.955e8 | m
Mjupiter = 1.8987e27 | kg
Rjupiter = 71492.0 | km
Mearth = 5.9722e24 | kg
Rearth = 6371.0088 | km  # IUGG mean radius
