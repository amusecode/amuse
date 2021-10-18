/* Standalone c program to use viscdisk to run the Pringle (1981) ring
   problem, with a complex equation of state */

#include <math.h>
#include <stdio.h>
#include "init.h"
#include "vader.h"
#include "vader_common.h"
#include <gsl/gsl_sf_bessel.h>

/* Physical constants */
#define KB              1.3806488e-16
#define MH              1.6737236e-24
#define MU              0.61
#define A               7.5657314e-15

/* Problem parameters */

/* Set this to a non-empty value to restart from a checkpoint */
#define CHECKPOINT      "ringrad_00001.vader"

/* Grid parmeters */
#define NR              4096
#define RMIN            1.5e10
#define RMAX            1.5e12
#define LINEAR          1

/* EOS parameters */
#define EOS_FUNC        1
#define GAMMA           1.66666666667
#define DELTA           0.0

/* Viscosity parameters */
#define ALPHA_FUNC      1
#define ALPHA           0.0

/* Inner boundary condition */
#define IBC_PRES_TYPE   FIXED_TORQUE
#define IBC_ENTH_TYPE   FIXED_ENTHALPY_VALUE
#define IBC_FUNC        1
#define IBC_PRES_VAL    0.0
#define IBC_ENTH_VAL    0.0

/* Outer boundary condition */
#define OBC_PRES_TYPE   FIXED_TORQUE
#define OBC_ENTH_TYPE   FIXED_ENTHALPY_VALUE
#define OBC_FUNC        1
#define OBC_PRES_VAL    0.0
#define OBC_ENTH_VAL    0.0

/* Source functions */
#define MASS_SRC_FUNC   0
#define MASS_SRC_VAL    0.0
#define INT_EN_SRC_FUNC 0
#define INT_EN_SRC_VAL  0.0

/* Control parameters */
#define ERR_TOL         1.0e-6
#define DT_TOL          0.1
#define MAX_ITER        40
#define INTERP_ORDER    1
#define MAX_DT_INCREASE 1.5
#define DT_MIN          1.0e-20
#define MAX_STEP        100000
#define USE_BE          0
#define PRE_TIMESTEP    0
#define POST_TIMESTEP   0
#define VERBOSITY       1

/* Output parameters */
#define NSAVE            65
#define OUTFILE         "ringrad.out"
#define NUSEROUT        0
#define CHECKNAME       "ringrad"
#define USERREADCHK     false
#define USERWRITECHK    false

/* Time parameters */
#define DT_START        -1
#define END_TIME        0.128

/* Problem-specific parameters*/
#define MSTAR           (1.99e33*1.4)
#define RING_MASS       1.99e27
#define INIT_TEMP       1.0e4
#define COL_RATIO       1.0e10
#define RING_LOC        7.5e11
#define NU              5.93e8
#define FZ0             7.5e9
#define GAMMA_GAS       (5.0/3.0)

/* Main */
int main () {

  grid *grd;
  double col0, t0, colAnalyt, x, tau;
  double *col, *pres, *eInt;
  double *colOut, *presOut, *mBndOut, *eBndOut, *eIntOut, *tOut;
  double param[5];
  double tSave[NSAVE];
  bool writeCheckpoint[NSAVE];
  char *checkpoint = CHECKPOINT;
  int i, j, idx;
  FILE *fp;
  unsigned long nStep, nIter, nFail, nOut;

  /* Allocate grid and workspace */
  grd = gridInitKeplerian(NR, RMIN, RMAX, MSTAR, LINEAR);

  /* Allocate data */
  col = (double *) calloc(NR, sizeof(double));
  pres = (double *) calloc(NR, sizeof(double));
  eInt = (double *) calloc(NR, sizeof(double));

  /* Initialize column density and pressure arrays, and associated variables */
  col0 = RING_MASS / (M_PI*SQR(RING_LOC));
  t0 = SQR(RING_LOC)/(12.0*NU);
  for (i=0, idx=-1; i<NR; i++)
    if ((grd->r_h[i] < RING_LOC) && (RING_LOC <= grd->r_h[i+1])) idx = i;
  for (i=0; i<NR; i++)
    col[i] = RING_MASS / grd->area[idx] / COL_RATIO;
  col[idx] = RING_MASS / grd->area[idx];
  for (i=0; i<NR; i++) {
    pres[i] = col[i]*INIT_TEMP*KB/(MU*MH) + (1./3.)*FZ0*A*SQR(SQR(INIT_TEMP));
    eInt[i] = col[i]*INIT_TEMP*KB/(MU*MH*(GAMMA_GAS-1)) + 
      FZ0*A*SQR(SQR(INIT_TEMP));
  }
  param[0] = NU;
  param[1] = RING_LOC;
  param[2] = RING_MASS;
  param[3] = RING_MASS / grd->area[idx] / COL_RATIO;
  param[4] = INIT_TEMP*KB/(MU*MH);

  /* Set output times */
  for (i=0; i<NSAVE; i++) {
    tSave[i] = t0 * i*END_TIME/((float) NSAVE-1);
    writeCheckpoint[i] = true;
  }

  /* Run the simulation */
  vader(
	/* Restart file name; if this is empty or NULL, calculation is
	   run from the start */
	checkpoint,
	/* Start and end time */
	0.0, END_TIME*t0,
	/* Equation of state parameters */
	EOS_FUNC, GAMMA, DELTA,
	/* Viscosity parameters */
	ALPHA_FUNC, ALPHA,
	/* Inner boundary condition */
	IBC_PRES_TYPE, IBC_ENTH_TYPE,
	IBC_FUNC, IBC_PRES_VAL, IBC_ENTH_VAL,
	/* Outer boundary condition */
	OBC_PRES_TYPE, OBC_ENTH_TYPE,
	OBC_FUNC, OBC_PRES_VAL, OBC_ENTH_VAL,
	/* Source functions */
	MASS_SRC_FUNC, MASS_SRC_VAL,
	INT_EN_SRC_FUNC, INT_EN_SRC_VAL,
	/* Control and method parameters */
	DT_START, DT_MIN, DT_TOL, ERR_TOL, MAX_DT_INCREASE,
	MAX_ITER, INTERP_ORDER, MAX_STEP, USE_BE,
	PRE_TIMESTEP, POST_TIMESTEP, VERBOSITY,
	/* Output control parameters */
	NSAVE, tSave, NUSEROUT, NULL,
	writeCheckpoint, CHECKNAME,
	USERREADCHK, USERWRITECHK,
	/* Computational grid */
	&grd,
	/* Initial conditions, holders for final conditions */
	&col, &pres, &eInt,
	/* User-defined extra parameters */
	&param,
	/* Diagnostic outputs */
	&nStep, &nIter, &nFail, &nOut,
	/* Simulation outputs */
	&tOut, &colOut, &presOut, &eIntOut, &mBndOut, &eBndOut,
	NULL, NULL, NULL
	);
	
  /* Print diagnostic output */
  printf("Total iterations = %ld, failed convergences = %ld\n",
	 nIter, nFail);

  /* Write output to file */
  if (!(fp=fopen(OUTFILE, "w"))) {
    fprintf(stderr, "Unable to open output file %s!\n", OUTFILE);
    exit(1);
  }
  fprintf(fp, "time/ts        x              col            colExact       pres           eInt\n");
  for (i=0; i<nOut; i++) {
    for (j=0; j<NR; j++) {
      x = grd->r_g[j+1]/RING_LOC;
      tau = tOut[i]/t0;
      colAnalyt = col0/ (pow(x, 0.25)*tau) * exp(-SQR(x-1.0)/tau) *
	gsl_sf_bessel_Inu_scaled(0.25, 2*x/tau);
      if (!isfinite(colAnalyt)) colAnalyt=0.0;
      if (colAnalyt < RING_MASS / grd->area[idx] / COL_RATIO)
	colAnalyt = RING_MASS / grd->area[idx] / COL_RATIO;
      fprintf(fp, "%e   %e   %e   %e   %e   %e\n", 
	      tau, x,
	      colOut[grd->nr*i+j], 
	      colAnalyt,
	      presOut[grd->nr*i+j],
	      eIntOut[grd->nr*i+j]);
    }
  }
  fclose(fp);

  /* Return successful exit */
  return(0);
}



