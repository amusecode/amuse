#ifndef GALACTICS_MAIN_COMMON_H
#define GALACTICS_MAIN_COMMON_H

/* Common global variables used by all executables */

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

#endif

