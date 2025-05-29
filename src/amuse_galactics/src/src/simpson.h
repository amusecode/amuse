#ifndef GALACTICS_SIMPSON_H
#define GALACTICS_SIMPSON_H

/* Integrate the function 'fcn' using Simpson's rule */
/*  x0, xn - limits of integrand
 		 n - number of divisions used in Simpson's rule */

float simpson(float (*fcn)(float), float x0, float xn, int n);

#endif

