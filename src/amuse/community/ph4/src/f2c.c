
/* The ugly code nstab.c was created from nstab.f (as provided by
 * Mardling, without any modification), using f2c.  This file contains
 * a few minimal functions from the f2c library needed to make C/C++
 * compilation work.
 */

#include <math.h>
double pow_dd(double *x, double *y) {return pow(*x, *y);}
float pow_ri(float *r, int *i) {return pow(*r, *i);}
double pow_di(double *x, int *i) {return pow(*x, *i);}
