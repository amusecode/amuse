
/* The ugly code nstab.c was created from nstab.f (as provided by
 * Mardling, without any modification), using f2c.  This file contains
 * minimal definitions from the real f2c.h to make C/C++ compilation
 * work.
 */

/* #include "f2c.h" */

typedef float		real;
typedef double		doublereal;
typedef long int	integer;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
