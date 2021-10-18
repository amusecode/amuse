/* Standalone c program to use viscdisk to run the Krumholz & Burkert
   (2010, ApJ, 724, 895) test problem. */

#include <math.h>
#include <stdio.h>
#include "init.h"
#include "driver.h"
#include "vader_common.h"

/* Physical constants */
#define G                 6.67384e-8   /* 2010 CODATA value in CGS units */

/* Problem parameters */
/* Grid parameters */
#define NR                512
#define RMIN              3.09e20
#define RMAX              3.09e22
#define LINEAR            0

/* EOS parameters */
#define EOS_FUNC          0
#define GAMMA             (5./3.)
#define DELTA             0.0

/* Viscosity parameters */
#define ALPHA_FUNC        1
#define ALPHA_VAL         0.0

/* Inner boundary condition parameters */
#define IBC_PRES_TYPE   FIXED_TORQUE
#define IBC_ENTH_TYPE   FIXED_ENTHALPY_VALUE
#define IBC_FUNC        1
#define IBC_PRES_VAL    0.0
#define IBC_ENTH_VAL    0.0

/* Outer boundary condition parameters */
#define OBC_PRES_TYPE   FIXED_TORQUE
#define OBC_ENTH_TYPE   FIXED_ENTHALPY_VALUE
#define OBC_FUNC        0.0
#define OBC_PRES_VAL    -6.30e25    /* 1 Msun/yr */

/* Source functions */
#define MASS_SRC_FUNC   0
#define MASS_SRC_VAL    0.0
#define INT_EN_SRC_FUNC 1
#define INT_EN_SRC_VAL  0.0

/* Control parameters */
#define ERR_TOL         1.0e-6
#define DT_TOL          0.1
#define MAX_ITER        40
#define INTERP_ORDER    2
#define MAX_DT_INCREASE 1.5
#define DT_MIN          1.0e-15
#define MAX_STEP        1
#define USE_BE          0
#define VERBOSITY       3

/* Output parameters */
#define NOUT            41
#define OUTFILE         "output/gidisk.out"

/* Time parameters */
#define N_ORBIT         20.0

/* Problem-specific parameters */
#define VPHI            2.2e7
#define MDOT            (-OBC_PRES_VAL)
#define ETA             1.5
#define T_Q             1.0
#define INIT_COL        2.0
#define INIT_VDISP      2.0
#define OBC_VDISP       2.0
#define DT_INIT         1.0e-5

/* Main */
int main() {

  grid *grd;
  wksp *w;
  double chi, s, torb, h_steady, obc_enth_val, dt_start, col1;
  int i, j;
  unsigned long nStep, nIter, nFail;
  double param[3];
  double col[NR], pres[NR], colSteady[NR], presSteady[NR];
  double *tOut, *colOut, *presOut, *QOut;
  double *mBndOut, *eBndOut, *eSrcOut;
  FILE *fp;

  /* Allocate grid and workspace */
  grd = gridInitFlat(NR, RMIN, RMAX, VPHI, LINEAR);
  w = wkspAlloc(NR);

  /* Compute characteristic values */
  chi = G*MDOT/(VPHI*VPHI*VPHI);
  s = pow(chi/ETA, 1./3.)/sqrt(2);
  torb = 2.0*M_PI*RMAX/VPHI;

  /* Set enthalpies correctly */
  h_steady = GAMMA / (GAMMA-1) * SQR(s*VPHI);
  obc_enth_val = h_steady * SQR(OBC_VDISP) / INIT_COL;

  /* Set starting timestep */
  dt_start = DT_INIT * torb;

  /* Set parameters */
  param[0] = ETA;
  param[1] = chi;
  param[2] = T_Q;

  /* Set initial conditions and calculate steady state */
  col1 = VPHI*VPHI * pow(chi/ETA, 1./3.) / (M_PI*G*RMAX);
  for (i=0; i<NR; i++) {
    colSteady[i] = col1 * RMAX/grd->r_g[i+1];
    presSteady[i] = colSteady[i] * SQR(s*VPHI);
    col[i] = colSteady[i] * INIT_COL;
    pres[i] = presSteady[i] * INIT_COL * SQR(INIT_VDISP);
  }

  /* Allocate memory to store results */
  if (!(tOut = calloc(NOUT, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(colOut = calloc(NOUT*NR, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(presOut = calloc(NOUT*NR, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(QOut = calloc(NOUT*NR, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(eSrcOut = calloc(NOUT*NR, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(mBndOut = calloc(NOUT*2, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(eBndOut = calloc(NOUT*2, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }

  /* Set output times */
  for (i=0; i<NOUT; i++) tOut[i] = i*N_ORBIT*torb/(NOUT-1);

  /* Run the simulation */
  printf("obc enth = %e\n", obc_enth_val);
  driver(0, N_ORBIT*torb, dt_start, grd, 
	 col, pres, NULL,
	 EOS_FUNC, GAMMA, DELTA,
	 ALPHA_FUNC, ALPHA_VAL,
	 IBC_PRES_TYPE, IBC_ENTH_TYPE,
	 IBC_FUNC, IBC_PRES_VAL, IBC_ENTH_VAL,
	 OBC_PRES_TYPE, OBC_ENTH_TYPE,
	 OBC_FUNC, OBC_PRES_VAL, obc_enth_val,
	 MASS_SRC_FUNC, MASS_SRC_VAL,
	 INT_EN_SRC_FUNC, INT_EN_SRC_VAL,
	 ERR_TOL, DT_TOL, MAX_ITER,
	 INTERP_ORDER, MAX_DT_INCREASE, DT_MIN,
	 MAX_STEP, USE_BE, VERBOSITY,
	 w, &param,
	 NOUT, tOut, colOut, presOut, NULL,
	 mBndOut, eBndOut, NULL, eSrcOut,
	 &nStep, &nIter, &nFail);

  /* Print diagnostic output */
  printf("Total iterations = %ld, failed convergences = %ld\n",
	 nIter, nFail);

  /* Compute Q from output */
  for (j=0; j<NOUT; j++) {
    for (i=0; i<grd->nr; i++) {
      QOut[i+j*grd->nr] = sqrt(2.0*(grd->beta_g[i+1]+1)) *
	grd->vphi_g[i+1]/grd->r_g[i+1] *
	sqrt(presOut[i+j*grd->nr]/colOut[i+j*grd->nr]) /
	(M_PI*G*colOut[i+j*grd->nr]);
    }
  }

  /* Write output to file */
  if (!(fp=fopen(OUTFILE, "w"))) {
    fprintf(stderr, "Unable to open output file %s!\n", OUTFILE);
    exit(1);
  }
  fprintf(fp, "time/tOrb      r/R            col/col1       pres/pres1                Q   colSteady/col1 presSteady/pres1\n");
  for (j=0; j<NOUT; j++) {
    for (i=0; i<grd->nr; i++) {
      fprintf(fp, "%e   %e   %e   %e   %e   %e   %e\n",
	      tOut[j]/torb, grd->r_g[i+1]/RMAX,
	      colOut[i+grd->nr*j]/col1, presOut[i+grd->nr*j]/(col1*SQR(VPHI)),
	      QOut[i+grd->nr+j], colSteady[i]/col1,
	      presSteady[i]/(col1*SQR(VPHI)));
    }
  }

  /* Return successful exit */
  return(0);
}
