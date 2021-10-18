/**************************************************************************/
/* Common definitions for the c vader routines                            */
/**************************************************************************/


/**************************************************************************/
/* General note on naming conventions: arrays without subscripts have     */
/* indices that range from 0 to nr-1, with no ghost zones, and are cell-  */
/* centered. Arrays with _h are edge centered and have indices that run   */
/* from 0 to nr. Arrays with _g are cell-centered with one ghost zone,    */
/* and have indices that go from 0 to nr+1; indices 0 and nr are the      */
/* ghost zones.                                                           */
/**************************************************************************/

#ifndef _vader_common_h_
#define _vader_common_h_

#include <gsl/gsl_vector.h>
#include <float.h>
#include <stdbool.h>
#include <time.h>

/* Slope limit parameter */
#define SLOPELIMIT  0.1

/* Descriptor for a grid */
typedef struct {
  unsigned long nr;            /* Number of real cells */
  bool linear;                 /* Is this a linear or logarithmic grid? */
  double *r_g, *r_h;           /* Cell center, edge locations */
  double *dr_g;                /* Cell sizes / log sizes */
  double *area;                /* Area of a zone */
  double *vphi_g, *vphi_h;     /* Rotation curve */
  double *beta_g, *beta_h;     /* Logarithmic index of rotation curve */
  double *psiEff_g, *psiEff_h; /* Effective gravitational potential */
  double *g_h;                 /* Factor appearing in derivatives */
} grid;

/* Workspace for calculations */
typedef struct {
  double *pres_g, *presNew_g, *colNew, *colTmp;
  double *alpha_g, *hint_g, *hintL_g, *hintR_g;
  double *ppmwksp_g;
  double *fmLast_h, *fmNew_h;
  double *ftLast_h, *feLast_h;
  double *massSrcLast, *massSrcNew, *intEnSrc;
  double *eIntTmp, *eIntNew;
  double *gammaLast, *deltaLast, *gammaNew, *deltaNew;
  double *mSrc, *eSrc;
  gsl_vector *ud_g, *ld_g, *diag_g, *rhs_g, *presTmp_g;
#if AA_M > 0
  double *colHist, *presHist, *eIntHist;
  double *colResid, *presResid, *eIntResid;
  gsl_vector *constraint;
#endif
} wksp;

/* Pressure boundary condition types */
typedef enum { FIXED_MASS_FLUX, FIXED_TORQUE_FLUX, FIXED_TORQUE } 
  pres_bc_type;

/* Enthalpy boundary condition types */
typedef enum { FIXED_ENTHALPY_VALUE, FIXED_ENTHALPY_GRADIENT } 
  enth_bc_type;

/* IO status indicators */
typedef enum { GOOD_IO, IO_ERROR, ALLOCATION_ERROR } iostatus;

/* Startup status indicators */
typedef enum { GOOD_START, RESTART_ERROR, MEMORY_ERROR, FIRST_DT_ERROR }
  setup_status;

/* Simulation status indicators */
typedef enum { RUNNING, NORMAL_EXIT, ZENO_ERROR, TOO_MANY_STEPS }
  status;

/* Macros used various places in code */
#define SQR(x) ((x)*(x))
#define LARGE DBL_MAX
#define SMALL DBL_MIN

#endif 
/* end _vader_common_h_ */
