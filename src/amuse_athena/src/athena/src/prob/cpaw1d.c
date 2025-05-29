#include "copyright.h"
/*============================================================================*/
/*! \file cpaw1d.c
 *  \brief Problem generator for 1-D circularly polarized Alfven wave (CPAW) 
 *  test.
 *
 * PURPOSE: Problem generator for 1-D circularly polarized Alfven wave (CPAW)
 *   test.  Only works in 1D (wavevector in x).  Tests in 2D and 3D are 
 *   initialized with different functions.
 *
 *   Can be used for either standing (problem/v_par=1.0) or travelling
 *   (problem/v_par=0.0) waves.
 *
 * USERWORK_AFTER_LOOP function computes L1 error norm in solution by comparing
 *   to initial conditions.  Problem must be evolved for an integer number of
 *   wave periods for this to work.
 *
 * REFERENCE: G. Toth,  "The div(B)=0 constraint in shock capturing MHD codes",
 *   JCP, 161, 605 (2000)						      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : The cpaw1d test only works for mhd.
#endif

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
static ConsS *RootSoln=NULL;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js;
  int k, ks = pGrid->ks;
  ConsS *Soln;
  int dir;      /* right (1) or left (2) polarization */
  Real x1,x2,x3,cs,sn,b_par,b_perp,fac,lambda,k_par,v_par,v_perp,den,pres;
#ifdef RESISTIVITY
  Real v_A, kva, omega_h, omega_l, omega_r;
#endif

  if ((Soln = (ConsS*)malloc(((ie-is+1)+2*nghost)*sizeof(ConsS))) == NULL)
    ath_error("[problem] Error initializing Soln array");

  if (pDomain->Level == 0) {
    if ((RootSoln = (ConsS*)malloc(((ie-is+1)+2*nghost)*sizeof(ConsS)))
       == NULL) ath_error("[problem] Error initializing RootSoln array");
  }

  if (pGrid->Nx[1] > 1 || pGrid->Nx[2] > 1) {
    ath_error("[cpaw1d] grid must be 1D");
  }

/* Put one wavelength on the root Domain, and initialize k_parallel */

  lambda = pDomain->RootMaxX[0] - pDomain->RootMinX[0]; 
  k_par = 2.0*PI/lambda;

  b_par = par_getd("problem","b_par");
  den = 1.0;
  b_perp = par_getd("problem","b_perp");
  v_perp = b_perp/sqrt((double)den);

  dir    = par_geti_def("problem","dir",1); /* right(1)/left(2) polarization */
  if (dir == 1) /* right polarization */
    fac = 1.0;
  else          /* left polarization */
    fac = -1.0;

#ifdef RESISTIVITY
  Q_Hall = par_getd("problem","Q_H");
  d_ind  = 0.0;
  v_A    = b_par/sqrt((double)den);
  if (Q_Hall > 0.0)
  {
    kva     = k_par*v_A;
    omega_h = 1.0/Q_Hall; /* characteristic frequency for Hall */

    omega_r = 0.5*SQR(kva)/omega_h*(sqrt(1.0+SQR(2.0*omega_h/kva)) + 1.0);
    omega_l = 0.5*SQR(kva)/omega_h*(sqrt(1.0+SQR(2.0*omega_h/kva)) - 1.0);

    if (dir == 1) /* right polarization (whistler wave) */
      v_perp = v_perp * kva / omega_r;
    else          /* left polarization */
      v_perp = v_perp * kva / omega_l;
  }
#endif

/* The gas pressure and parallel velocity are free parameters. */
  pres = par_getd("problem","pres");
  v_par = par_getd("problem","v_par");

/* Setup circularily polarized AW solution  */

  for (i=is; i<=ie+1; i++) {
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);

    sn = sin(k_par*x1);
    cs = cos(k_par*x1);

    Soln[i].d = den;

    Soln[i].M1 = den*v_par;
    Soln[i].M2 = -fac*den*v_perp*sn;
    Soln[i].M3 = -den*v_perp*cs;
 
    Soln[i].B1c = b_par;
    Soln[i].B2c = fac*b_perp*sn;
    Soln[i].B3c = b_perp*cs;
#ifndef ISOTHERMAL
    Soln[i].E = pres/Gamma_1 
      + 0.5*den*(v_par*v_par + v_perp*v_perp)
      + 0.5*(b_par*b_par + b_perp*b_perp);
#endif
  }

/* set code variables to solution */

  for (i=is; i<=ie+1; i++) {
    pGrid->U[ks][js][i].d  = Soln[i].d;
    pGrid->U[ks][js][i].M1 = Soln[i].M1;
    pGrid->U[ks][js][i].M2 = Soln[i].M2;
    pGrid->U[ks][js][i].M3 = Soln[i].M3;

    pGrid->U[ks][js][i].B1c = pGrid->B1i[ks][js][i] = Soln[i].B1c;
    pGrid->U[ks][js][i].B2c = pGrid->B2i[ks][js][i] = Soln[i].B2c;
    pGrid->U[ks][js][i].B3c = pGrid->B3i[ks][js][i] = Soln[i].B3c;
#ifndef ISOTHERMAL
    pGrid->U[ks][js][i].E = Soln[i].E;
#endif
  }

/* save solution on root grid */

  if (pDomain->Level == 0) {
    for (i=is; i<=ie+1; i++) {
      RootSoln[i].d  = Soln[i].d;
      RootSoln[i].M1 = Soln[i].M1;
      RootSoln[i].M2 = Soln[i].M2;
      RootSoln[i].M3 = Soln[i].M3;

      RootSoln[i].B1c = Soln[i].B1c;
      RootSoln[i].B2c = Soln[i].B2c;
      RootSoln[i].B3c = Soln[i].B3c;
#ifndef ISOTHERMAL
      RootSoln[i].E = Soln[i].E;
#endif
    }
  }

  free(Soln);
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

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
/*! \fn void get_eta_user(GridS *pG, int i, int j, int k,
 *                           Real *eta_O, Real *eta_H, Real *eta_A)
 *  \brief Get user defined resistivity. */
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{
  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif


void Userwork_in_loop(MeshS *pM)
{
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in CPAW,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(MeshS *pM)
{
  GridS *pGrid;
  int i=0,is,ie,js,ks,Nx1;
  Real rms_error=0.0;
  ConsS error;
  FILE *fp;
  char *fname;

  error.d = 0.0;
  error.M1 = 0.0;
  error.M2 = 0.0;
  error.M3 = 0.0;
  error.B1c = 0.0;
  error.B2c = 0.0;
  error.B3c = 0.0;
#ifndef ISOTHERMAL
  error.E = 0.0;
#endif /* ISOTHERMAL */

/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js;
  ks = pGrid->ks;
  Nx1 = (ie-is+1);

/* compute L1 error in each variable, and rms total error */

  for (i=is; i<=ie; i++) {
    error.d   += fabs(pGrid->U[ks][js][i].d   - RootSoln[i].d  );
    error.M1  += fabs(pGrid->U[ks][js][i].M1  - RootSoln[i].M1 );
    error.M2  += fabs(pGrid->U[ks][js][i].M2  - RootSoln[i].M2 );
    error.M3  += fabs(pGrid->U[ks][js][i].M3  - RootSoln[i].M3 );
    error.B1c += fabs(pGrid->U[ks][js][i].B1c - RootSoln[i].B1c);
    error.B2c += fabs(pGrid->U[ks][js][i].B2c - RootSoln[i].B2c);
    error.B3c += fabs(pGrid->U[ks][js][i].B3c - RootSoln[i].B3c);
#ifndef ISOTHERMAL
    error.E   += fabs(pGrid->U[ks][js][i].E   - RootSoln[i].E  );
#endif /* ISOTHERMAL */
  }

/* Compute RMS error over all variables */

  rms_error = SQR(error.d) + SQR(error.M1) + SQR(error.M2) + SQR(error.M3);
  rms_error += SQR(error.B1c) + SQR(error.B2c) + SQR(error.B3c);
#ifndef ISOTHERMAL
  rms_error += SQR(error.E);
#endif /* ISOTHERMAL */
  rms_error = sqrt(rms_error)/(double)Nx1;

/* Print error to file "cpaw1d-errors.dat" */

  fname = "cpaw1d-errors.dat";
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_after_loop]: Unable to open file.\n");
      return;
    }
/* write out some header information */
    fprintf(fp,"# Nx1 RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
    fprintf(fp,"  B1c  B2c  B3c");
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %e  %e  %e  %e  %e",Nx1,rms_error,
          (error.d/(double)Nx1),
          (error.M1/(double)Nx1),
          (error.M2/(double)Nx1),
          (error.M3/(double)Nx1));
#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(error.E/(double)Nx1));
#endif /* ISOTHERMAL */
  fprintf(fp,"  %e  %e  %e",
          (error.B1c/(double)Nx1),
          (error.B2c/(double)Nx1),
          (error.B3c/(double)Nx1));
  fprintf(fp,"\n");

  fclose(fp);

  return;
}
