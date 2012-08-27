#include "copyright.h"
/*============================================================================*/
/*! \file pgflow.c
 *  \brief Problem generator for steady planar gravitational flow in a simple
 *   1D gravitational field
 *
 * PURPOSE: Problem generator for steady planar gravitational flow in a simple
 *   1D gravitational field: g = grav*cos(k_par*x) with periodic boundary
 *   conditions.  The 1D flow can be initialized in a 3D (x1,x2,x3) domain
 *   using the following transformation rules:
 *   -  x =  x1*cos(alpha) + x2*sin(alpha)
 *   -  y = -x1*sin(alpha) + x2*cos(alpha)
 *   -  z = x3
 *
 *   This problem is a good test of the source terms in a static grav potential.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - rtbis()          - finds roots via bisection
 * - grav_pot()       - gravitational potential
 * - Bfunc()          - computes Bernoilli function
 * - expr_drho()      - computes difference d-d0
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef ADIABATIC
#error The Bernoulli solver assumes an adiabatic eos.
#endif

#ifdef MHD
#error The Bernoulli solver assumes a hydrodynamic fluid.
#endif

static Real grav, psi;
static Real H, S, Phi; /* Bernoulli const., specific entropy, mass flux */
static Real sin_a, cos_a; /* sin and cos of alpha */
static Real lambda, k_par; /* Wavelength, 2.0*PI/wavelength */
static Real E0; /* The total initial energy (including the potential)
		   averaged over the computational grid. */
static int root; /* 0 -> super-sonic root, otherwise -> sub-sonic root */
static Real ***d0=NULL;  /* initial density, used by expr_drho */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * rtbis()          - finds roots via bisection
 * grav_pot()       - gravitational potential
 * Bfunc()          - computes Bernoilli function
 * expr_drho()      - computes difference d-d0
 *============================================================================*/

static int rtbis(double (*pfun)(double), const double x1, const double x2,
		 const double xacc, const int imax, double *prt);
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static double Bfunc(double rho);
static Real expr_drho(const Grid *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i, is = pGrid->is, ie = pGrid->ie, nx1;
  int j, js = pGrid->js, je = pGrid->je, nx2;
  int k, ks = pGrid->ks, ke = pGrid->ke, nx3;
  Real x1,x2,x3;
  Real den,v_par,pres;
  Real rho_p,rho_e,rho_s;
  Real dt;
  double rho;
  double angle; /* Angle the wave direction makes with the x1-direction */
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  if ((d0 = (Real***)calloc_3d_array(nx3,nx2,nx1,sizeof(Real))) == NULL)
    ath_error("[pgflow]: Error allocating memory\n");

/* An angle = 0.0 is a wave aligned with the x1-direction. */
  angle = par_getd("problem","angle");

/* If the grid is 1-D we override the angle variable */
  if(pGrid->Nx2 <= 1) angle = 0.0;
  if(pGrid->Nx1 <= 1) angle = 90.0;

/* Compute the sin and cos of the angle and the wavelength. */
  if (angle == 0.0) {
    sin_a = 0.0;
    cos_a = 1.0;
    lambda = pGrid->Nx1*pGrid->dx1; /* Put one wavelength in the grid */
  }
  else if (angle == 90.0) {
    sin_a = 1.0;
    cos_a = 0.0;
    lambda = pGrid->Nx2*pGrid->dx2; /* Put one wavelength in the grid */
  }
  else {
/* We put 1 wavelength in each direction.  Hence the wavelength
 *     lambda = pGrid->Nx1*pGrid->dx1*cos_a;
 *     AND  lambda = pGrid->Nx2*pGrid->dx2*sin_a;
 *     are both satisfied. */
    if((pGrid->Nx1*pGrid->dx1) == (pGrid->Nx2*pGrid->dx2)){
      cos_a = sin_a = sqrt(0.5);
    }
    else{
      angle = atan((double)(pGrid->Nx1*pGrid->dx1)/(pGrid->Nx2*pGrid->dx2));
      sin_a = sin(angle);
      cos_a = cos(angle);
    }
/* Use the larger angle to determine the wavelength */
    if (cos_a >= sin_a) {
      lambda = pGrid->Nx1*pGrid->dx1*cos_a;
    } else {
      lambda = pGrid->Nx2*pGrid->dx2*sin_a;
    }
  }

/* Initialize k_parallel */
  k_par = 2.0*PI/lambda;

  E0 = 0.0; /* Initialize the total energy */
  grav = par_getd("problem","grav");
  root = par_geti("problem","root");

  den = par_getd("problem","den");
  pres = par_getd("problem","pres");
  v_par = par_getd("problem","v_par");

/* Set up the constants of motion */
  Phi = den*v_par; /* The mass flux */
  S   = pres/pow((double)den,(double)Gamma); /* specific entropy */
/* Calculate the Bernoulli constant at x1 = 0.0 */
  H   = 0.5*v_par*v_par + Gamma*pres/(Gamma_1*den);

  rho_e = pow((double)(Phi*Phi/(Gamma*S)),(double)(1.0/(Gamma+1.0)));

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Calculate the gravitational potential */
      psi = -grav*sin((double)k_par*(x1*cos_a + x2*sin_a))/k_par;

      if(H <= psi)
	ath_error("[problem]: H < Psi -- No solution exists\n");

      if(Bfunc((double)rho_e) < 0.0)
	ath_error("[problem]: Bfunc(rho_e) < 0.0 -- No solution exists\n");

      if(root){ /* Choose the heavy (subsonic) root */
/* The root is bounded: rho_e <= rho < rho_s */
	rho_s = pow((double)(Gamma_1*(H-psi)/(Gamma*S)),(double)(1.0/Gamma_1));

	if(rtbis(Bfunc,rho_e,rho_s,1.0e-12*rho_e,100,&rho)){
	  exit(1);
	}
      } else { /* Choose the light (supersonic) root */
/* The root is bounded: rho_p < rho <= rho_e */
	rho_p = fabs((double)Phi)/sqrt((double)(2.0*(H-psi)));

	if(rtbis(Bfunc,rho_p,rho_e,1.0e-12*rho_e,100,&rho)){
	  exit(1);
	}
      }

      d0[k][j][i] = rho;
      pGrid->U[k][j][i].d  = rho;
      pGrid->U[k][j][i].M1 = Phi*cos_a;
      pGrid->U[k][j][i].M2 = Phi*sin_a;
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = S*pow(rho,(double)Gamma)/Gamma_1
        + 0.5*Phi*Phi/rho;

      E0 += pGrid->U[k][j][i].E + rho*psi;
    }
  }}

/* Average over the domain */
  E0 /= (Real)((ie - is + 1)*(je - js + 1)*(ke - ks + 1));

/* Enroll the gravitational potential function */
  StaticGravPot = grav_pot;

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
  if(strcmp(expr,"drho")==0) return expr_drho;
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn int rtbis(double (*pfun)(double), const double x1, const double x2,
 *		  const double xacc, const int imax, double *prt)
 *  \brief  Using bisection, find the root of a function "pfun" known to lie
 * between x1 and x2 to an accuracy <= "xacc", but do not exceed
 * "imax" bisections.  
 *
 * The root is returned through the pointer "prt".
 * RETURN VALUE: 0 on Success, 1 on Error
 * ASSUMPTIONS: This routine assumes that a first order zero lies
 * between x1 and x2, i.e. the function is continuous between x1 and
 * x2 and in a sufficiently small neighborhood of the root, the
 * function is monotonic. Written by T. A. Gardiner -- Sept. 24, 2003 
 */

int rtbis(double (*pfun)(double), const double x1, const double x2,
	  const double xacc, const int imax, double *prt)
{
  int i;
  double xn, xm, xp, dx;
  double fn = (*pfun)(x1), fm, fp = (*pfun)(x2);

  if(fn < 0.0 && fp > 0.0){
    xn = x1;
    xp = x2;
  }
  else if(fn > 0.0 && fp < 0.0){
    xn = x2;
    xp = x1;

    fm = fn;
    fn = fp;
    fp = fm;
  }
  else if(fn == 0.0){ *prt = x1;  return 0; }
  else if(fp == 0.0){ *prt = x2;  return 0; }
  else{
    ath_perr(-1,"[rtbis]: Root must be bracketed for bisection\n");
    ath_perr(-1,"[rtbis]: x1=%g, x2=%g, F(x1)=%g, F(x2)=%g\n",
	     x1,x2,fn,fp);
    return 1; /* Error */
  }

  dx = xp - xn;

  for(i=0; i<imax; i++){ /* Bisection loop */
    dx *= 0.5;
    xm = xn + dx;
    fm = (*pfun)(xm);

    if(fm < 0.0){ xn = xm;  fn = fm; }
    else        { xp = xm;  fp = fm; }

    if(fabs(dx) < xacc || fm == 0.0){
      *prt = fp < -fn ? xp : xn; /* Return our best value for the root */
      return 0;
    }
  }

  ath_perr(-1,"[rtbis]: Too many bisections\n");

  *prt = fp < -fn ? xp : xn; /* Return our best value for the root */

  return 1;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential
 */

static Real grav_pot(const Real x1, const Real x2, const Real x3)
{
  return -grav*sin((double)k_par*(x1*cos_a + x2*sin_a))/k_par;
}

/*----------------------------------------------------------------------------*/
/*! \file static double Bfunc(double rho)
 *  \brief Computes Bernoulli function
 */

static double Bfunc(double rho){
  return (H - psi - 0.5*Phi*Phi/(rho*rho) 
	  - (Gamma*S/Gamma_1)*pow((double)rho,(double)Gamma_1));
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_drho(const Grid *pG,const int i,const int j,
 *			      const int k)
 *  \brief Computes d-d0 (density - initial value)
 */

static Real expr_drho(const Grid *pG, const int i, const int j, const int k)
{
 return (pG->U[k][j][i].d - d0[k][j][i]);
}
