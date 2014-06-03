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


/* Disk Constants */
float mdisk, rdisk, zdisk;
float rd, zd;
float vsigR, vsigp, vsigz, vphimax;
float diskmass, diskedge;

/* Bulge Constants */
float rho1=100, sigb=0.3, sigb2, psicutb= -1.0, fbulgeconst;
float bulgemass, bulgeedge;
int bulgeflag;


/* Halo Constants */
float A, B, C, D1, D2, D3, q2, v0, v02;
float halomass, haloedge;
int haloflag;

