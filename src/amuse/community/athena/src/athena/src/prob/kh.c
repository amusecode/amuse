#include "copyright.h"
/*============================================================================*/
/*! \file kh.c
 *  \brief Problem generator for KH instability. 
 *
 * PURPOSE: Problem generator for KH instability.  Sets up two versions:
 * - iprob=1: slip surface with random perturbations
 * - iprob=2: tanh profile at interface, with single-mode perturbation	      */
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
 * ran2()
 *============================================================================*/

static double ran2(long int *idum);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real amp,drat,vflow,b0,a,sigma,x1,x2,x3;
  long int iseed = -1;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read problem parameters */

  iprob = par_geti("problem","iprob");
  vflow = par_getd("problem","vflow");
  drat = par_getd("problem","drat");
  amp = par_getd("problem","amp");
#ifdef MHD
  b0  = par_getd("problem","b0");
#endif

/* iprob=1.  Two uniform streams moving at +/- vflow, random perturbations */

  if (iprob == 1) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          pGrid->U[k][j][i].d = 1.0;
          pGrid->U[k][j][i].M1 = vflow + amp*(ran2(&iseed) - 0.5);
          pGrid->U[k][j][i].M2 = amp*(ran2(&iseed) - 0.5);
          pGrid->U[k][j][i].M3 = 0.0;
          if (fabs(x2) < 0.25) {
  	    pGrid->U[k][j][i].d = drat;
            pGrid->U[k][j][i].M1 = -drat*(vflow + amp*(ran2(&iseed) - 0.5));
            pGrid->U[k][j][i].M2 = drat*amp*(ran2(&iseed) - 0.5);
          }
/* Pressure scaled to give a sound speed of 1 with gamma=1.4 */
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E = 2.5/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* BAROTROPIC */
#ifdef MHD
          pGrid->B1i[k][j][i] = b0;
          pGrid->U[k][j][i].B1c = b0;
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif /* BAROTROPIC */
#endif /* MHD */
        }
#ifdef MHD
      pGrid->B1i[k][j][ie+1] = b0;
#endif
      }
    }
  }

/* iprob=2.  Test suggested by E. Zweibel, based on Ryu & Jones.
 * Two uniform density flows with single mode perturbation
 */

  if (iprob == 2) {
    a = 0.05;
    sigma = 0.2;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          pGrid->U[k][j][i].d = 1.0;
          pGrid->U[k][j][i].M1 = vflow*tanh(x2/a);
          pGrid->U[k][j][i].M2 = amp*sin(2.0*PI*x1)*exp(-(x2*x2)/(sigma*sigma));
          pGrid->U[k][j][i].M3 = 0.0;
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E = 1.0/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* BAROTROPIC */
#ifdef MHD
          pGrid->B1i[k][j][i] = b0;
          pGrid->U[k][j][i].B1c = b0;
#ifndef BAROTROPIC
          pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif /* BAROTROPIC */
#endif /* MHD */
/* Use passive scalar to keep track of the fluids, since densities are same */
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = 0.0;
          if (x2 > 0) pGrid->U[k][j][i].s[0] = 1.0;
#endif
        }
#ifdef MHD
      pGrid->B1i[k][j][ie+1] = b0;
#endif
      }
    }
  }

/* With viscosity and/or resistivity, read eta_R and nu_V */

#ifdef OHMIC
  eta_Ohm = par_getd("problem","eta");
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif
#ifdef BRAGINSKII
  nu_V = par_getd("problem","nu");
#endif

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

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
#ifdef OHMIC
  eta_Ohm = par_getd("problem","eta");
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif
#ifdef BRAGINSKII
  nu_V = par_getd("problem","nu");
#endif
  return;
}

#if (NSCALARS > 0)
/*! \fn static Real color(const GridS *pG, const int i, const int j,const int k)
 *  \brief Returns first passively advected scalar s[0] */
static Real color(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"color")==0) return color;
#endif
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*------------------------------------------------------------------------------
 */

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
 *  \brief  Extracted from the Numerical Recipes in C (version 2) code. Modified
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
double ran2(long int *idum){
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
