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
  float chalo;
  float v0;
  float a;
  float nnn;
  float v0bulge;
  float abulge;
  float psi0;
  float haloconst;
  float bulgeconst;
  float mdisk;
  float rdisk;
  float zdisk;
  float outdisk;
  float drtrunc;
  float potcor;
} galaxy;

galaxy gparam;

typedef struct {
  float psic;
  float psid;
} energytable;

energytable cparam;

typedef struct {
  float fcut_halo;
  float fcut_bulge;
} dfcutoff;

dfcutoff bparam;

phase *r;


/* Disk Constants */
float mdisk, rdisk, zdisk, outdisk, drtrunc;
float rd, zd;
float vsigR, vsigp, vsigz, vphimax;
float diskmass, diskedge;

/* Bulge Constants */
float nnn, v0bulge, abulge, bulgeconst;
float bulgemass, bulgeedge;
float psic, psid;
float fcut_halo;

/* Halo Constants */

float chalo, v0, a, haloconst, v02, psi0;
float halomass, haloedge;
float fcut_bulge;

/* Blackhole Constants */
float bhmass, bhsoft;
