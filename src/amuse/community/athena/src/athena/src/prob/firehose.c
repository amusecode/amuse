#include "copyright.h"
/*============================================================================*/
/*! \file firehose.c
 *  \brief Problem generator for firehose test of Braginskii viscosity.
 *
 * PURPOSE: Problem generator for firehose test of Braginskii viscosity.
 * - iprob=1: vortex test suggested by Steve Cowley, developed by Greg Hammett
 * - iprob=2: shear test from Appendix C.3.4 in Prateek Sharma's thesis.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - ran2()    - random number generator from NR
 * - pbc_ix1() - sets BCs on L-x1 (left edge) of grid used in shear test
 * - pbc_ox1() - sets BCs on R-x1 (right edge) of grid used in shear test */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()    - random number generator from NR
 * pbc_ix1() - sets BCs on L-x1 (left edge) of grid used in shear test
 * pbc_ox1() - sets BCs on R-x1 (right edge) of grid used in shear test
 *============================================================================*/

static double ran2(long int *idum);
static void pbc_ix1(Grid *pGrid);
static void pbc_ox1(Grid *pGrid);

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i,is,ie,j,js,je,ks,nx1,nx2,iprob;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real d0,p0,B0,v0,x1,x2,x3,x1min,x1max,Lx,k0,press,amp,kp;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  nx1 = (ie-is)+1;
  nx2 = (je-js)+1;
  if ((nx1 == 1) || (nx2 == 1)) {
    ath_error("[firehose]: This problem can only be run in 2D\n");
  }
  d0 = 1.0;
  p0 = 1.0;
  B0 = par_getd("problem","B0");
  v0 = par_getd("problem","v0");
#ifdef BRAGINSKII
  nu_V = par_getd("problem","nu");
#endif
  amp = par_getd("problem","amp");
  kp =  par_getd("problem","kp");
  iprob =  par_getd("problem","iprob");

/* Initialize wavenumber */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx = x1max - x1min;
  k0 = 2.0*PI/Lx;

/* iprob=1: Cowley/Hammett vortex test ---------------------------------------*/
/* Initialize density, momentum, face-centered fields */

  if (iprob == 1) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);

      pGrid->U[ks][j][i].d = d0;
      pGrid->U[ks][j][i].M1 = -d0*v0*cos(k0*x2)*sin(k0*x1);
      pGrid->U[ks][j][i].M2 = d0*v0*sin(k0*x2)*cos(k0*x1)+amp*cos(kp*x1);
      pGrid->U[ks][j][i].M3 = 0.0;
#ifdef MHD
      pGrid->B1i[ks][j][i] = B0;
      pGrid->B2i[ks][j][i] = 0.0;
      pGrid->B3i[ks][j][i] = 0.0;
      pGrid->U[ks][j][i].B1c = B0;
      pGrid->U[ks][j][i].B2c = 0.0;
      pGrid->U[ks][j][i].B3c = 0.0;
#endif /* MHD */
      press = p0 + 0.5*d0*v0*v0*(cos(k0*x1)*cos(k0*x1) + cos(k0*x2)*cos(k0*x2));
      pGrid->U[ks][j][i].E = press/Gamma_1
#ifdef MHD
        + 0.5*(SQR(pGrid->U[ks][j][i].B1c) + SQR(pGrid->U[ks][j][i].B2c)
             + SQR(pGrid->U[ks][j][i].B3c))
#endif /* MHD */
        + 0.5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M2)
              + SQR(pGrid->U[ks][j][i].M3))/pGrid->U[ks][j][i].d;
    }
  }
#ifdef MHD
/* boundary conditions on interface B */
  for (j=js; j<=je; j++) {
    pGrid->B1i[ks][j][ie+1] = B0;
  }
  for (i=is; i<=ie; i++) {
    pGrid->B2i[ks][je+1][i] = 0.0;
  }
#endif /* MHD */
  }

/* iprob=2: Sharma shear test ----------------------------------------------- */
/* Initialize density, momentum, face-centered fields */

  if (iprob == 2) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);

      pGrid->U[ks][j][i].d = d0;
      pGrid->U[ks][j][i].M1 = d0*amp*(ran2(&iseed) - 0.5);
      pGrid->U[ks][j][i].M2 = d0*(amp*(ran2(&iseed) - 0.5) - 0.015*x1);
      pGrid->U[ks][j][i].M3 = 0.0;
#ifdef MHD
      pGrid->B1i[ks][j][i] = B0;
      pGrid->B2i[ks][j][i] = B0;
      pGrid->B3i[ks][j][i] = 0.0;
      pGrid->U[ks][j][i].B1c = B0;
      pGrid->U[ks][j][i].B2c = B0;
      pGrid->U[ks][j][i].B3c = 0.0;
#endif /* MHD */
      pGrid->U[ks][j][i].E = 0.1/Gamma_1
#ifdef MHD
        + 0.5*(SQR(pGrid->U[ks][j][i].B1c) + SQR(pGrid->U[ks][j][i].B2c)
             + SQR(pGrid->U[ks][j][i].B3c))
#endif /* MHD */
        + 0.5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M2)
              + SQR(pGrid->U[ks][j][i].M3))/pGrid->U[ks][j][i].d;
    }
  }
#ifdef MHD
/* boundary conditions on interface B */
  for (j=js; j<=je; j++) {
    pGrid->B1i[ks][j][ie+1] = B0;
  }
  for (i=is; i<=ie; i++) {
    pGrid->B2i[ks][je+1][i] = B0;
  }
#endif /* MHD */
/* Enroll special BCs for shearing sheet */
  set_bvals_mhd_fun(left_x1,  pbc_ix1);
  set_bvals_mhd_fun(right_x1, pbc_ox1);
  }

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3,
                                          Real *w1, Real *w2, Real *w3)
{
  return;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

/*----------------------------------------------------------------------------*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief Extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 * 
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

/*----------------------------------------------------------------------------*/
/*! \fn static void pbc_ix1(Grid *pGrid)
 *  \brief Special PERIODIC boundary conditions, Inner x1 boundary
 */

static void pbc_ix1(Grid *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
  Real xmin,xmax,Lx;
  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][is-i].M2 =
          pGrid->U[k][j][is-i].M2 + pGrid->U[k][j][is-i].d*0.015*Lx;
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][is-i] = pGrid->W[k][j][ie-(i-1)];
#endif
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pbc_ox1(Grid *pGrid)
 *  \brief Special PERIODIC boundary conditions, Outer x1 boundary
 */

static void pbc_ox1(Grid *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
  Real xmin,xmax,Lx;
  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][ie+i].M2 =
          pGrid->U[k][j][ie+i].M2 - pGrid->U[k][j][ie+i].d*0.015*Lx;
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][ie+i] = pGrid->W[k][j][is+(i-1)];
#endif
      }
    }
  }

#ifdef MHD
/* Note that i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}
