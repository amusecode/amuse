#ifndef GALACTICS_MAIN_H
#define GALACTICS_MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define GCON 6.349363592e-2 /* (2.0*pi)^(-3/2) */

/* Generate a distribution which follows an oblate logarithmic potential */

typedef struct {
	float x, y, z, vx, vy, vz;
}  phase;

typedef struct {
	float A;
	float B;
	float C;
	float v0;
	float q;
	float psi0;
	float psicut;
	float rho1;
	float sigbulge;
	float mdisk;
	float rdisk;
	float zdisk;
	float routdisk;
	float drtrunc;
	float potcor;
	int idiskflag;
	int ibulgeflag;
	int ihaloflag;
} galaxy;

phase *r;
galaxy gparam;

#endif

