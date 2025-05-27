#include "copyright.h"
/*============================================================================*/
/*! \file cshock1d.c
 *  \brief Problem generator for 1-D standing C-type shock test.
 *
 * PURPOSE: Problem generator for 1-D standing C-type shock test. Only works for
 *   grid-aligned shock, the shock is along the x1 direction.
 *
 * REFERENCE: Mac Low, M-M et al., "Incorporation of Ambipolar Diffusion into
 *   the Zeus Magnetohydrodynamics Code", 1995, ApJ, 442, 726		      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error: The cshock1d test only works for mhd.
#endif

#ifndef ISOTHERMAL
#error: The cshock1d test only works for isothermal equation of state.
#endif

#ifndef RESISTIVITY
#error: The cshock1d test only works with RESISTIVITY.
#endif

/*--------------- Private Functions --------------------*/
/* RK4 integrator for the semi-analytical solution */

Real RK4(Real D, Real A, Real M, Real theta, Real h);

Real Dprime(Real D, Real A, Real M, Real theta);

/* Solution in the root domain */
static ConsS *RootSoln=NULL;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  ConsS *Soln;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int i0, nx1, N_L;
  Real x1min,x1max,Lx,x1,x2,x3,xs,xe,h,x01,x02;
  Real Mach,Alfv,theta,d0,v0,vA,B0,Bx0,By0;
  Real Ls,Ns,D01,D02,myD;

  nx1 = (ie-is)+1 + 2*nghost;

/* allocate memory for solution on this level.  If this is root level
 * also allocate memory for RootSoln */
  if ((Soln = (ConsS*)calloc_1d_array(nx1,sizeof(ConsS)))==NULL)
    ath_error("[problem]: Error allocating memory for Soln\n");
  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS*)calloc_1d_array(nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }

/* Model parameters */
  Mach = par_getd("problem","Mach"); /* Mach number */
  Alfv = par_getd("problem","Alfv"); /* Alfven Mach number */

  theta = par_getd("problem","theta"); /* magnetic field obliquety (deg) */
  theta = theta * PI / 180.0;

/* upstream quantities (in the shock frame) */
  d0 = 1.0;
  v0 = Mach * Iso_csound;

  vA = (Mach/Alfv) * Iso_csound;
  B0 = sqrt(SQR(vA)*d0);

  Bx0 = B0 * cos(theta);
  By0 = B0 * sin(theta);

/* root domain info */
  x1min = pDomain->RootMinX[0];
  x1max = pDomain->RootMaxX[0];
  Lx = x1max - x1min;

/* Ambipolar diffusion coefficient 
 * ambipolar diffusion length scale is FIXED to be 1 in code unit! */

  Q_AD = 1.0/vA;            /* N.B. AD diffusion coefficient is set this way! */
  d_ind= 0.0;               /* constant ion density */

/* info for semi-analytic calculation */
  Ls  = par_getd_def("problem","Ls",20.0);  /* Estimated thickness of the
                                             * C-shock (in unit of L_A) */
  Ns  = par_getd_def("problem","Ns",5e3);   /* number of cells in Ls in the
                                             * semi-analytic calculation */
  if (Ls > Lx)
    ath_error("Domain size in x1 is shorter than the C-Shosck thickness!\n");

  xs = x1min + 0.5*(Lx - Ls);
  xe = xs + Ls;
  h  = (xe - xs)/Ns;

/* Find the 1D solution on this grid.
 * Semi-analytic solution is only applied in the middle of the root domain.
 * Have to initialize different segments separately. */

  /* Upstream region */
  i = is;
  cc_pos(pGrid,i,js,ks,&x1,&x2,&x3); 
  while (x1 < xs){

    Soln[i].d  = d0;
    Soln[i].M1 = d0 * v0;
    Soln[i].M2 = 0.0;
    Soln[i].M3 = 0.0;
    Soln[i].B1c= Bx0;
    Soln[i].B2c= By0;
    Soln[i].B3c= 0.0;

    i += 1;
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
  }

  /* C-shock region */
  x01 = xs;		x02 = x01+h;
  D01 = d0+1.0e-6;
  while (x02 <= xe)
  {
    D02 = RK4(D01, Alfv, Mach, theta, h);
    if ((x1 >= x01) && (x1 < x02))
    {
      myD = (D01*(x02-x1)+D02*(x1-x01))/h;
      Soln[i].d  = myD;
      Soln[i].B1c= Bx0;
      Soln[i].B2c= sqrt(SQR(By0)+2.0*SQR(Alfv*B0)*
                       (myD-1)*(1.0/myD-SQR(1.0/Mach)));
      Soln[i].B3c= 0.0;
      Soln[i].M1 = d0 * v0;
      Soln[i].M2 = myD* SQR(vA)/v0*cos(theta)*(Soln[i].B2c/B0-sin(theta));
      Soln[i].M3 = 0.0; 

      i += 1;
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
    }
    x01 = x02;
    x02 += h;
    D01 = D02;
  }
  i0 = i-1;

  /* Downstream region */
  while (i <= ie+1)
  {
    Soln[i] = Soln[i0];
    i += 1;
  }

  /* backup the root domain solution */
  if (pDomain->Level == 0){
    for (i=is; i<=ie+1; i++){
      RootSoln[i].d  = Soln[i].d;
      RootSoln[i].M1 = Soln[i].M1;
      RootSoln[i].M2 = Soln[i].M2;
      RootSoln[i].M3 = Soln[i].M3;
      RootSoln[i].B1c= Soln[i].B1c;
      RootSoln[i].B2c= Soln[i].B2c;
      RootSoln[i].B3c= Soln[i].B3c;
    }
  }

  /* set initial conditions */
  ie = ie+1;
  if (pGrid->Nx[1] > 1) je = je+1;
  if (pGrid->Nx[2] > 1) ke = ke+1;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {

    pGrid->U[k][j][i].d  = Soln[i].d;
    pGrid->U[k][j][i].M1 = d0 * v0;
    pGrid->U[k][j][i].M2 = Soln[i].M2;
    pGrid->U[k][j][i].M3 = 0.0;

    pGrid->U[k][j][i].B1c = pGrid->B1i[k][j][i] = Bx0;
    pGrid->U[k][j][i].B2c = pGrid->B2i[k][j][i] = Soln[i].B2c;
    pGrid->U[k][j][i].B3c = pGrid->B3i[k][j][i] = 0.0;
  }}}

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
}

/*----------------------------------------------------------------------------*/
/* Semi-analytical C-shock solution */

/*! \fn Real RK4(Real D, Real A, Real M, Real theta, Real h)
 *  \brief 4th order Runge-Kutta integrator */
Real RK4(Real D, Real A, Real M, Real theta, Real h)
{
  Real k1, k2, k3, k4;

  k1 = Dprime(D,A,M,theta);
  k2 = Dprime(D+0.5*h*k1,A,M,theta);
  k3 = Dprime(D+0.5*h*k2,A,M,theta);
  k4 = Dprime(D+h*k3,A,M,theta);

  return D + h/6.0*(k1+2.0*k2+2.0*k3+k4);

}

/*! \fn Real Dprime(Real D, Real A, Real M, Real theta)
 *  \brief dD/dt */
Real Dprime(Real D, Real A, Real M, Real theta)
{
  Real sintheta, costheta, sintheta2, costheta2;
  Real M21, b, b2;

  sintheta = sin(theta);
  costheta = cos(theta);
  sintheta2 = SQR(sintheta);
  costheta2 = SQR(costheta);
  M21 = 1.0/SQR(M);

  b2 = sintheta2 + 2*SQR(A)*(D-1.0)*(1.0/D-M21);
  b  = sqrt(b2);

  return b/A*(b-D*((b-sintheta)/SQR(A)*costheta2+sintheta))/(b2+costheta2)/(1/SQR(D)-M21);
}

