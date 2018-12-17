/* Standard headers */
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
/* User-defined header */
#include "physical_constants.h"

namespace physical_constants {

// Physical constants
const double Ggrav = 6.67408e-8;  // Gravitational constants [cm^{3}/g/s^{2}]
const double kBoltz = 1.3806503e-16; // Boltzmann constant [cm^{2} g/s^{2}/K]

// Mass units
const double Msolar = 1.989e33; // Solar mass [g]
const double dalton = 1.660538921e-24; // Unified atomic mass unit [g]
const double Mhydrogen = 1.007825 * dalton; // Mass of hydrogen atom [g]

// Length units
const double km = 1.0e5; // kilo-meters [cm]
const double AU = 1.49597870e13; // Astronomical unit [cm]
const double pc = 3.0857e18; // parsec [cm]
const double kpc = 1.0e3 * pc; // kilo-parsecs [cm]
const double Mpc = 1.0e6 * pc; // mega-parsecs [cm]
const double Gpc = 1.0e9 * pc; // giga-parsecs [cm]

// Time units
const double yr = 3.15576e7; // year [s]
const double kyr = 1.0e3 * yr; // kilo-years [s]
const double Myr = 1.0e6 * yr; // mega-years [s]
const double Gyr = 1.0e9 * yr; // giga-years [s]

}

