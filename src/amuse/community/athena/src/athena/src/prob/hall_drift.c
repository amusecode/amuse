#include "copyright.h"
/*==============================================================================
 * FILE: hall_drift.c
 *
 * PURPOSE: Problem generator for Hall drift mode test in 2D. The problem
 *   is initiated by a non-linear sinusoidal distribution of B_z field in the x1
 *   direction and a sinusoidal distribution of electron density in the x2
 *   direction. The system evolves as Burger's equation in x1, with shear in x2.
 *
 *   Configure --enable-resistivity
 *
 *   IMPORTANT NOTICE:
 *     To perform this test, the main integrator MUST be turned off!
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : The hall drift test only works for mhd.
#endif

#ifndef RESISTIVITY
#error : The hall drift test only works with resistivity.
#endif

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3,Nx1,Nx2;
  Real x1size,x2size;

  Real B0,dB,d0,dden,Pres; 
  Real k1, k2;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

/* NOTE: For parallel calculations Nx1 != nx1 and Nx2 != nx2 */

  Nx1 = pDomain->Nx[0];
  Nx2 = pDomain->Nx[1];
  if (Nx1 == 1 || Nx2 == 1) {
    ath_error("[problem]: this test only works with Nx1 & Nx2 > 1\n");
  }

/* Read field strength and fluctuations */
  B0   = par_getd_def("problem","B0",1);
  dB   = par_getd_def("problem","dB",0.1);
  d0   = 1.0;
  dden = par_getd_def("problem","drho",0.1);
  Pres = SQR(0.5*B0);

/* Set angle, wavelength */

  x1size = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  x2size = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  k1 = 2.0*PI/x1size;
  k2 = 2.0*PI/x2size;

/* Now initialize 2D solution */

#ifdef MHD
/* Set face-centered fields */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie+1; i++) {
    pGrid->B1i[k][j][i] = 0.0;
  }}}
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->B2i[k][j][i] = 0.0;
  }}}
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    Real x1,x2,x3;
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    pGrid->B3i[k][j][i] = B0 + dB*cos(k1*x1);
  }}}
  if (pGrid->Nx[2] > 1) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->B3i[ke+1][j][i] = pGrid->B3i[ke][j][i];
    }}
  }
#endif

/* Now set initial conditions to 2d wave solution */ 

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    Real x1,x2,x3;
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    pGrid->U[k][j][i].d = d0/(1.0 - dden*cos(k2*x2)/d0);

    pGrid->U[k][j][i].M1 = 0.0;
    pGrid->U[k][j][i].M2 = 0.0;
    pGrid->U[k][j][i].M3 = 0.0;
#ifdef MHD
    pGrid->U[k][j][i].B1c = 0.0;
    pGrid->U[k][j][i].B2c = 0.0;
    pGrid->U[k][j][i].B3c = B0 + dB*cos(k1*x1);
#endif /* MHD */

#ifndef ISOTHERMAL
    pGrid->U[k][j][i].E = Pres/Gamma_1 +
      0.5*pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c;
#endif /* ISOTHERMAL */

  }}}

#ifdef RESISTIVITY 
  eta_Ohm = 0.0;
  Q_AD    = 0.0;
  Q_Hall  = par_getd("problem","Q_H");
  d_ind   = 1.0;
#endif

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
#ifdef RESISTIVITY  
  eta_Ohm = 0.0;
  Q_AD    = 0.0;
  Q_Hall  = par_getd("problem","Q_H");
  d_ind   = 1.0;
#endif
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

void Userwork_after_loop(MeshS *pM)
{
}
