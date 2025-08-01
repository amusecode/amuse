#include "copyright.h"
/*============================================================================*/
/*! \file cpaw2d.c
 *  \brief Problem generator for 2-D circularly polarized Alfven wave (CPAW) 
 *  test.
 *
 * PURPOSE: Problem generator for 2-D circularly polarized Alfven wave (CPAW)
 *   test.  Works for any arbitrary wavevector in the x1-x2 plane.  The wave is
 *   defined with reference to a coordinate system (x,y,z) with transformation
 *   rules to the code coordinate system (x1,x2,x3)
 *   -  x =  x1*cos(alpha) + x2*sin(alpha)
 *   -  y = -x1*sin(alpha) + x2*cos(alpha)
 *   -  z = x3
 *   The magnetic field is given by:
 *   - B_x = b_par
 *   - B_y = b_perp*sin(k*x)
 *   - B_z = b_perp*cos(k*x)   where k = 2.0*PI/lambda
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
#error : The cpaw2d test only works for MHD.
#endif


/* Initial solution, shared with Userwork_after_loop to compute L1 error.
 * Vector potential A3 defined in prvate function below.  B_par, etc. must
 * be defined as globals to be used by A3() */
 
static ConsS **RootSoln=NULL;
static Real A3(const Real x1, const Real x2);
Real fac, sin_a, cos_a, b_par, b_perp;
Real k_par;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  ConsS **Soln;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int nx1, nx2;
  int dir;
  Real angle;    /* Angle the wave direction makes with the x1-direction */
  Real x1size,x2size,x1,x2,x3,cs,sn;
  Real v_par, v_perp, den, pres;
  Real lambda; /* Wavelength */
#ifdef RESISTIVITY
  Real v_A, kva, omega_h, omega_l, omega_r;
#endif
 
  nx1 = (ie - is + 1) + 2*nghost;
  nx2 = (je - js + 1) + 2*nghost;

  if (pGrid->Nx[1] == 1) {
    ath_error("[problem] Grid must be 2D");
  }

  if ((Soln = (ConsS**)calloc_2d_array(nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[problem]: Error allocating memory for Soln\n");

  if (pDomain->Level == 0){
    if ((RootSoln =(ConsS**)calloc_2d_array(nx2,nx1,sizeof(ConsS)))==NULL)
      ath_error("[problem]: Error allocating memory for RootSoln\n");
  }

/* An angle =  0.0 is a wave aligned with the x1-direction. */
/* An angle = 90.0 is a wave aligned with the x2-direction. */

  angle = par_getd("problem","angle");
  x1size = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  x2size = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

/* Compute the sin and cos of the angle and the wavelength. */
/* Put one wavelength in the grid */

  if (angle == 0.0) {
    sin_a = 0.0;
    cos_a = 1.0;
    lambda = x1size;
  }
  else if (angle == 90.0) {
    sin_a = 1.0;
    cos_a = 0.0;
    lambda = x2size;
  }
  else {

/* We put 1 wavelength in each direction.  Hence the wavelength
 *      lambda = (pDomain->RootMaxX[0] - pDomain->RootMinX[0])*cos_a
 *  AND lambda = (pDomain->RootMaxX[1] - pDomain->RootMinX[1])*sin_a;
 *  are both satisfied. */

    if(x1size == x2size){
      cos_a = sin_a = sqrt(0.5);
    }
    else{
      angle = atan((double)(x1size/x2size));
      sin_a = sin(angle);
      cos_a = cos(angle);
    }
/* Use the larger angle to determine the wavelength */
    if (cos_a >= sin_a) {
      lambda = x1size*cos_a;
    } else {
      lambda = x2size*sin_a;
    }
  }

/* Initialize k_parallel */

  k_par = 2.0*PI/lambda;
  b_par = par_getd("problem","b_par");
  den = 1.0;

  ath_pout(0,"va_parallel = %g\nlambda = %g\n",b_par/sqrt(den),lambda);

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
    omega_h = 1.0/Q_Hall;

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

/* Use the vector potential to initialize the interface magnetic fields
 * The iterface fields are located at the left grid cell face normal */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        cs = cos(k_par*(x1*cos_a + x2*sin_a));

        x1 -= 0.5*pGrid->dx1;
        x2 -= 0.5*pGrid->dx2;

        pGrid->B1i[k][j][i] = -(A3(x1,(x2+pGrid->dx2)) - A3(x1,x2))/pGrid->dx2;
        pGrid->B2i[k][j][i] =  (A3((x1+pGrid->dx1),x2) - A3(x1,x2))/pGrid->dx1;
        pGrid->B3i[k][j][i] = b_perp*cs;
      }
    }
  }
  if (pGrid->Nx[2] > 1) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        cs = cos(k_par*(x1*cos_a + x2*sin_a));
        pGrid->B3i[ke+1][j][i] = b_perp*cs;
      }
    }
  }

/* Now initialize the cell centered quantities */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

        sn = sin(k_par*(x1*cos_a + x2*sin_a)) * fac;
        cs = cos(k_par*(x1*cos_a + x2*sin_a));

        Soln[j][i].d  = den;
        Soln[j][i].M1 = den*(v_par*cos_a + v_perp*sn*sin_a);
        Soln[j][i].M2 = den*(v_par*sin_a - v_perp*sn*cos_a);
        Soln[j][i].M3 = -den*v_perp*cs;
        pGrid->U[k][j][i].d  = Soln[j][i].d;
        pGrid->U[k][j][i].M1 = Soln[j][i].M1;
        pGrid->U[k][j][i].M2 = Soln[j][i].M2;
        pGrid->U[k][j][i].M3 = Soln[j][i].M3;

        Soln[j][i].B1c = 0.5*(pGrid->B1i[k][j][i] + pGrid->B1i[k][j][i+1]);
        Soln[j][i].B2c = 0.5*(pGrid->B2i[k][j][i] + pGrid->B2i[k][j+1][i]);
        Soln[j][i].B3c = b_perp*cs;
        pGrid->U[k][j][i].B1c = Soln[j][i].B1c;
        pGrid->U[k][j][i].B2c = Soln[j][i].B2c;
        pGrid->U[k][j][i].B3c = Soln[j][i].B3c;

#ifndef ISOTHERMAL
        Soln[j][i].E = pres/Gamma_1 + 0.5*(SQR(pGrid->U[k][j][i].B1c) +
                 SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c) )
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) +
                 SQR(pGrid->U[k][j][i].M3) )/den;
        pGrid->U[k][j][i].E = Soln[j][i].E;
#endif
      }
    }
  }

/* save solution on root grid */

  if (pDomain->Level == 0) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      RootSoln[j][i].d  = Soln[j][i].d ;
      RootSoln[j][i].M1 = Soln[j][i].M1;
      RootSoln[j][i].M2 = Soln[j][i].M2;
      RootSoln[j][i].M3 = Soln[j][i].M3;
#ifndef ISOTHERMAL
      RootSoln[j][i].E  = Soln[j][i].E ;
#endif /* ISOTHERMAL */
#ifdef MHD
      RootSoln[j][i].B1c = Soln[j][i].B1c;
      RootSoln[j][i].B2c = Soln[j][i].B2c;
      RootSoln[j][i].B3c = Soln[j][i].B3c;
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        RootSoln[j][i].s[n] = Soln[j][i].s[n];
#endif
    }}
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
 * A3() - computes vector potential to initialize fields
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
  int i,j,is,ie,js,je,ks,Nx1,Nx2;
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
  if (pGrid == NULL) return;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  Nx1 = (ie-is+1);
  Nx2 = (je-js+1);

/* compute L1 error in each variable, and rms total error */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[ks][j][i].d   - RootSoln[j][i].d  );
      error.M1  += fabs(pGrid->U[ks][j][i].M1  - RootSoln[j][i].M1 );
      error.M2  += fabs(pGrid->U[ks][j][i].M2  - RootSoln[j][i].M2 );
      error.M3  += fabs(pGrid->U[ks][j][i].M3  - RootSoln[j][i].M3 );
      error.B1c += fabs(pGrid->U[ks][j][i].B1c - RootSoln[j][i].B1c);
      error.B2c += fabs(pGrid->U[ks][j][i].B2c - RootSoln[j][i].B2c);
      error.B3c += fabs(pGrid->U[ks][j][i].B3c - RootSoln[j][i].B3c);
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[ks][j][i].E   - RootSoln[j][i].E  );
#endif /* ISOTHERMAL */
    }
  }

/* Compute RMS error over all variables */

  rms_error = SQR(error.d) + SQR(error.M1) + SQR(error.M2) + SQR(error.M3);
  rms_error += SQR(error.B1c) + SQR(error.B2c) + SQR(error.B3c);
#ifndef ISOTHERMAL
  rms_error += SQR(error.E);
#endif /* ISOTHERMAL */
  rms_error = sqrt(rms_error)/(double)(Nx1*Nx2);

/* Print error to file "cpaw2d-errors.dat" */

  fname = "cpaw2d-errors.dat";
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
    fprintf(fp,"# Nx1   Nx2  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
    fprintf(fp,"  B1c  B2c  B3c");
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d %e  %e  %e  %e  %e",Nx1,Nx2,rms_error,
          (error.d/(double)(Nx1*Nx2)),
          (error.M1/(double)(Nx1*Nx2)),
          (error.M2/(double)(Nx1*Nx2)),
          (error.M3/(double)(Nx1*Nx2)));
#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(error.E/(double)(Nx1*Nx2)));
#endif /* ISOTHERMAL */
  fprintf(fp,"  %e  %e  %e",
          (error.B1c/(double)(Nx1*Nx2)),
          (error.B2c/(double)(Nx1*Nx2)),
          (error.B3c/(double)(Nx1*Nx2)));
  fprintf(fp,"\n");

  fclose(fp);

  return;
}

/*---------------------------------------------------------------------------*/
/*! \fn static Real A3(const Real x1, const Real x2)
 *  \brief Define a scalar potential A3 such that:
 *   - B_x = - $\partial A3 / \partial y$
 *   - B_y =   $\partial A3 / \partial x$
 *   Then A3 is given in the function below.  */

static Real A3(const Real x1, const Real x2)
{
  return b_par*(x1*sin_a - x2*cos_a) 
     - (fac*b_perp/k_par)*cos(k_par*(x1*cos_a + x2*sin_a));
}
