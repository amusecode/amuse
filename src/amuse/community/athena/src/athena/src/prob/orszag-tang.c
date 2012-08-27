#include "copyright.h"
/*============================================================================*/
/*! \file orszag-tang.c
 *  \brief Problem generator for Orszag-Tang vortex problem.
 *
 * REFERENCE: For example, see: G. Toth,  "The div(B)=0 constraint in shock
 *   capturing MHD codes", JCP, 161, 605 (2000)				      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : The problem generator orszag_tang.c only works for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,is,ie,j,js,je,ks,nx1,nx2;
  Real d0,p0,B0,v0,x1,x2,x3,**az;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  if (pGrid->Nx[2] > 1) {
    ath_error("[orszag-tang]: This problem can only be run in 2D\n");
  }
  if (pGrid->Nx[1] == 1) {
    ath_error("[orszag-tang]: This problem can only be run in 2D\n");
  }
  if ((az = (Real**)calloc_2d_array(nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[orszag-tang]: Error allocating memory for vector pot\n");
  }
  B0 = 1.0/sqrt(4.0*PI);
  d0 = 25.0/(36.0*PI);
  v0 = 1.0;
  p0 = 5.0/(12*PI);

/* Initialize vector potential */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
      x1 -= 0.5*pGrid->dx1;
      x2 -= 0.5*pGrid->dx2;
      az[j][i] = B0/(4.0*PI)*cos(4.0*PI*x1) + B0/(2.0*PI)*cos(2.0*PI*x2);
    }
  }

/* Initialize density, momentum, face-centered fields */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
/* Calculate the cell center positions */
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);

      pGrid->U[ks][j][i].d = d0;
      pGrid->U[ks][j][i].M1 = -d0*v0*sin(2.0*PI*x2);
      pGrid->U[ks][j][i].M2 =  d0*v0*sin(2.0*PI*x1);
      pGrid->U[ks][j][i].M3 = 0.0;
      pGrid->B1i[ks][j][i] = (az[j+1][i] - az[j][i])/pGrid->dx2;
      pGrid->B2i[ks][j][i] =-(az[j][i+1] - az[j][i])/pGrid->dx1;
    }
  }
/* boundary conditions on interface B */
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie+1; i++) {
    pGrid->B1i[ks][j][i] = (az[j+1][i] - az[j][i])/pGrid->dx2;
  }}
  for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->B2i[ks][j][i] =-(az[j][i+1] - az[j][i])/pGrid->dx1;
  }}

/* initialize total energy and cell-centered B */

  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->U[ks][j][i].B1c = 0.5*(pGrid->B1i[ks][j][i]+pGrid->B1i[ks][j][i+1]);
    pGrid->U[ks][j][i].B2c = 0.5*(pGrid->B2i[ks][j][i]+pGrid->B2i[ks][j+1][i]);
    pGrid->U[ks][j][i].B3c = 0.0;
#ifndef ISOTHERMAL
    pGrid->U[ks][j][i].E = p0/Gamma_1
        + 0.5*(SQR(pGrid->U[ks][j][i].B1c) + SQR(pGrid->U[ks][j][i].B2c)
             + SQR(pGrid->U[ks][j][i].B3c))
        + 0.5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M2)
              + SQR(pGrid->U[ks][j][i].M3))/pGrid->U[ks][j][i].d;
#endif
  }}

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

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
