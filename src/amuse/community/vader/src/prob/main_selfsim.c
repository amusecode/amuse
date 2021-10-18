/* Standalone c program to use viscdisk to run the Lynden-Bell &
   Pringle (1974, MNRAS, 168, 603) self-similar viscous disk test problem. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "init.h"
#include "vader.h"

/* Physical constants */
#define KB              1.3806488e-16
#define MH              1.6737236e-24
#define MU              2.33

/* Problem parameters */

/* Set this to a non-empty value to restart from a checkpoint */
#define CHECKPOINT      ""

/* Grid parmeters */
#define NR              512
#define RMIN            1.5e12
#define RMAX            3.0e14
#define LINEAR          0

/* EOS parameters */
#define EOS_FUNC        0
#define GAMMA           1.000001
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
#define INTERP_ORDER    2
#define MAX_DT_INCREASE 1.5
#define DT_MIN          1.0e-15
#define MAX_STEP        10000
#define USE_BE          0
#define PRE_TIMESTEP    0
#define POST_TIMESTEP   0
#ifndef TESTING_MODE
#   define VERBOSITY    1
#else
#   define VERBOSITY    0
#endif

/* Output parameters */
#define NSAVE           31
#define OUTFILE         "output/selfsim.out"
#define NUSEROUT        0
#define WRITECHECKPOINT NULL
#define CHECKNAME       NULL
#define USERREADCHK     false
#define USERWRITECHK    false

/* Time parameters */
#define DT_START        -1
#define END_TIME        4.0

/* Problem-specific parameters */
#define MSTAR           1.99e33
#define MDOT0           6.3e19
#define R0              1.5e13
#define NU0             2.37e13
#define INIT_TEMP       100.0


/* Main */
int main() {

  grid *grd;
  double col1, ts;
  double *col, *pres;
  double tSave[NSAVE], *colOut, *presOut, *tOut, *mBndOut, *eBndOut;
  double param[4];
  char *checkpoint = CHECKPOINT;
  FILE *fp;
  unsigned long nStep, nIter, nFail, nOut;
  int i, j;
#ifdef TESTING_MODE
  double *residSum;
  unsigned long *iterStep;
  double driverTime, advanceTime, nextIterTime, userTime;
#endif

  /* Allocate grid */
  grd = gridInitKeplerian(NR, RMIN, RMAX, MSTAR, LINEAR);

  /* Allocate data */
  col = (double *) calloc(NR, sizeof(double));
  pres = (double *) calloc(NR, sizeof(double));

#ifdef TESTING_MODE
  /* Allocate memory to hold output of tests */
  if (!(residSum = calloc(MAX_ITER, sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
  if (!(iterStep = calloc(MAX_STEP, sizeof(unsigned long)))) {
    fprintf(stderr, "Cannot allocate memory for output data!\n");
    exit(1);
  }
#endif
  
  /* Initialize column density and pressure arrays, and associated variables */
  col1 = MDOT0 / (3.0*M_PI*NU0);
  ts = R0*R0/(3*NU0);
  for (i=0; i<NR; i++) {
    col[i] = col1 * exp(-grd->r_g[i+1]/R0) / (grd->r_g[i+1]/R0);
    pres[i] = col[i] * INIT_TEMP*KB/(MU*MH);
  }
  param[0] = NU0;
  param[1] = R0;
  param[2] = MDOT0;
  param[3] = INIT_TEMP*KB/(MU*MH);

  /* Set save times */
  for (i=0; i<NSAVE; i++) 
    tSave[i] = ts * (1 + i*(END_TIME-1)/((float) NSAVE-1));

  /* Run the simulation */
  vader(
	/* Restart file name; if this is empty or NULL, calculation is
	   run from the start */
	checkpoint,
	/* Start and end time */
	ts, END_TIME*ts,
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
	WRITECHECKPOINT, CHECKNAME,
	USERREADCHK, USERWRITECHK,
	/* Computational grid */
	&grd,
	/* Initial conditions, holders for final conditions */
	&col, &pres, NULL,
	/* User-defined extra parameters */
	&param,
	/* Diagnostic outputs */
	&nStep, &nIter, &nFail, &nOut,
	/* Simulation outputs */
	&tOut, &colOut, &presOut, NULL, &mBndOut, &eBndOut,
	NULL, NULL, NULL
#ifdef TESTING_MODE
	, residSum, iterStep,
	&driverTime, &advanceTime, &nextIterTime, &userTime
#endif
	 );

  /* Print diagnostic output */
  printf("Total iterations = %ld, failed convergences = %ld\n",
	 nIter, nFail);
#ifdef TESTING_MODE
  printf("Timings (sec):\n");
  printf("   driver   : %f\n", driverTime);
  printf("   advance  : %f\n", advanceTime);
  printf("   nextIter : %f\n", nextIterTime);
  printf("   user     : %f\n", userTime);
#endif

  /* Write output to file */
  if (!(fp=fopen(OUTFILE, "w"))) {
    fprintf(stderr, "Unable to open output file %s!\n", OUTFILE);
    exit(1);
  }
  fprintf(fp, "time/ts        x              col/col0       colExact/col0  pres/pres0\n");
  for (i=0; i<NSAVE; i++) {
    for (j=0; j<NR; j++) {
      fprintf(fp, "%e   %e   %e   %e   %e\n", tOut[i]/ts, grd->r_g[j+1]/R0,
	      colOut[grd->nr*i+j]/col1, 
	      exp(-(grd->r_g[j+1]/R0)/(tOut[i]/ts)) /
	      (grd->r_g[j+1]/R0 * pow(tOut[i]/ts, 1.5)),
	      presOut[grd->nr*i+j]/(col1*INIT_TEMP*KB/(MU*MH)));
    }
  }
  fclose(fp);

  /* Free memory */
  free(col);
  free(pres);
  outputFree(&tOut, &colOut, &presOut, NULL,
	     &mBndOut, &eBndOut, NULL, NULL, NULL);
  gridFree(grd);
#ifdef TESTING_MODE
  free(residSum);
  free(iterStep);
#endif
  
  /* Return successful exit */
  return(0);
}
