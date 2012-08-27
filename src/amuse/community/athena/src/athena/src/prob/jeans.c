#include "copyright.h"
/*============================================================================*/
/*! \file jeans.c
 *  \brief Problem generator for simple self-gravity test.  
 *
 * PURPOSE: Problem generator for simple self-gravity test.  
 *
 * SELF_GRAVITY must be defined; use --with-gravity=fft 
 *               and for fft method, --enable-fft   
 *
 * B-field (when present) lies along direction perpendicular to wavevector.
 * Wavevector is along chosen direction (1, 2, or 3). */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = (pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n;
  int kmax,jmax;
  Real amp,njeans;
  int kdir;
  Real d0,p0,u0,v0,w0,b0;
  Real x1,x2,x3,sinkx,coskx,xl;
  Real lambda,omega,omega2,cs,va,kwave;
#ifdef MHD
  Real beta;
#endif /* MHD */
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

#ifndef SELF_GRAVITY
    ath_error("[jeans]: this test only works with SELF_GRAVITY configured");
#endif /* SELF_GRAVITY */

/* Read parameters for initial conditions */
  amp = par_getd("problem","amp");
#ifdef MHD
  beta = par_getd("problem","beta");
#endif
/* njeans = lambda/lambda_J */
  njeans = par_getd("problem","njeans");

/* direction =1, 2, or 3 */
  kdir =par_geti("problem","kdir");
/* wavelength = size of Domain along k direction*/
  switch(kdir){
    case 1:
    lambda = pDomain->Nx[0]*pGrid->dx1;
    break;
  case 2:
    lambda = pDomain->Nx[1]*pGrid->dx2;
    break;
  case 3:
    lambda = pDomain->Nx[2]*pGrid->dx3;
  }
/* background density =1 */
  d0 = 1.0;
/* background thermal pressure =1 */
  p0 = 1.0;
/* background velocities =0 */
  u0 = 0.0;
  v0 = 0.0;
  w0 = 0.0;
#ifdef MHD
  b0 = sqrt(2.*p0/beta);
  va=b0/sqrt(d0);
#else
  va = 0.0;
  b0 = 0.0;
#endif

#ifdef SELF_GRAVITY
/* Set gravity constant*/
  four_pi_G = (4.0*Gamma*p0)*(PI*PI*njeans*njeans)/(d0*d0*lambda*lambda);
  grav_mean_rho = d0;
#endif /* SELF_GRAVITY */
/* define useful parameters */
#ifndef ISOTHERMAL
  cs=sqrt(Gamma*p0/d0);
#else
    cs=Iso_csound;
#endif
    kwave=2.0*PI/lambda;
/* dispersion relation */
    omega2=kwave*kwave*(cs*cs + va*va - four_pi_G*d0/(kwave*kwave));
    omega=sqrt(fabs(omega2));

printf("4piG=%e, lambda=%e, period=%e\n",four_pi_G, lambda, (2.0*PI/omega));

 kmax= (pGrid->Nx[2]>1)? ke+1 : ke;
 jmax= (pGrid->Nx[1]>1)? je+1 : je;

/* Now initialize solution */
  for (k=ks; k<=kmax; k++) {
  for (j=js; j<=jmax; j++) {
  for (i=is; i<=ie+1; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
  switch(kdir){
    case 1:
    xl=x1-0.5*pGrid->dx1;
    sinkx = sin(x1*kwave);
    coskx = cos(x1*kwave);
    break;
  case 2:
    xl=x2-0.5*pGrid->dx2;
    sinkx = sin(x2*kwave);
    coskx = cos(x2*kwave);
    break;
  case 3:
    xl=x3-0.5*pGrid->dx3;
    sinkx = sin(x3*kwave);
    coskx = cos(x3*kwave);
  }
    pGrid->U[k][j][i].d = d0*(1.0 + amp*sinkx);
#ifndef ISOTHERMAL
    pGrid->U[k][j][i].E = (p0/Gamma_1)*(1.0+Gamma*amp*sinkx) + 
                b0*b0*(0.5+amp*sinkx);
#endif /* ISOTHERMAL */
  switch(kdir){
    case 1:
    pGrid->U[k][j][i].M1 = (omega2<0.0)? d0*(omega/kwave)*amp*coskx: 0.0;
    pGrid->U[k][j][i].M2 = 0.0;
    pGrid->U[k][j][i].M3 = 0.0;
    break;
  case 2:
    pGrid->U[k][j][i].M1 = 0.0;
    pGrid->U[k][j][i].M2 = (omega2<0.0)? d0*(omega/kwave)*amp*coskx: 0.0;
    pGrid->U[k][j][i].M3 = 0.0;
    break;
  case 3:
    pGrid->U[k][j][i].M1 = 0.0;
    pGrid->U[k][j][i].M2 = 0.0;
    pGrid->U[k][j][i].M3 = (omega2<0.0)? d0*(omega/kwave)*amp*coskx: 0.0;
  }
#ifdef MHD
  switch(kdir){
    case 1:
    pGrid->B1i[k][j][i] = 0.0;
    pGrid->B2i[k][j][i] = b0*(1.0+amp*sin(kwave*xl));
    pGrid->B3i[k][j][i] = 0.0;
    pGrid->U[k][j][i].B1c = 0.0;
    pGrid->U[k][j][i].B2c = b0*(1.0+amp*sinkx);
    pGrid->U[k][j][i].B3c = 0.0;
    break;
    case 2:
    pGrid->B1i[k][j][i] = b0*(1.0+amp*sin(kwave*xl));
    pGrid->B2i[k][j][i] = 0.0;
    pGrid->B3i[k][j][i] = 0.0;
    pGrid->U[k][j][i].B1c = b0*(1.0+amp*sinkx);
    pGrid->U[k][j][i].B2c = 0.0;
    pGrid->U[k][j][i].B3c = 0.0;
    break;
    case 3:
    pGrid->B1i[k][j][i] = b0*(1.0+amp*sin(kwave*xl));
    pGrid->B2i[k][j][i] = 0.0;
    pGrid->B3i[k][j][i] = 0.0;
    pGrid->U[k][j][i].B1c = b0*(1.0+amp*sinkx);
    pGrid->U[k][j][i].B2c = 0.0;
    pGrid->U[k][j][i].B3c = 0.0;
  }
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pGrid->U[k][j][i].s[n] = 0.0;
#endif
  }}}

  return;
}

/*=============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * color()   - returns first passively advected scalar s[0]
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

#if (NSCALARS > 0)
/*! \fn static Real color(const GridS *pG, const int i, const int j,const int k)
 *  \brief returns first passively advected scalar s[0] */
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
}


void Userwork_after_loop(MeshS *pM)
{
}
