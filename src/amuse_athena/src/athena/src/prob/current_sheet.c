#include "copyright.h"
/*============================================================================*/
/*! \file current_sheet.c
 *  \brief Problem generator for current sheet test. 
 *
 * PURPOSE: Problem generator for current sheet test.  This version only allows
 *   current sheet in X-Y plane, with Bz=0.  
 *
 * REFERENCE: */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef HYDRO
#error : The current sheet test only works for MHD.
#endif /* HYDRO */

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  Prim1DS W;
  Cons1DS U1d;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real x1,x2,x3;
  Real Bx,uflow,beta;

/* setup uniform ambient medium with current sheet */

  uflow = par_getd("problem","uflow");
  beta  = par_getd("problem","beta");
  W.d = 1.0;
  W.P = beta;
  W.Vy = 0.0;
  W.Vz = 0.0;
  Bx = 0.0;
  W.By = 1.0;
  W.Bz = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        W.Vx = uflow*cos(PI*x2);
        U1d = Prim1D_to_Cons1D(&(W), &Bx);
        pGrid->U[k][j][i].d  = U1d.d;
        pGrid->U[k][j][i].M1 = U1d.Mx; 
        pGrid->U[k][j][i].M2 = U1d.My;
        pGrid->U[k][j][i].M3 = U1d.Mz;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E = U1d.E;
#endif
        pGrid->U[k][j][i].B1c = Bx;
        pGrid->U[k][j][i].B2c = U1d.By;
        pGrid->U[k][j][i].B3c = U1d.Bz;
        pGrid->B1i[k][j][i] = Bx;
        pGrid->B2i[k][j][i] = U1d.By;
        pGrid->B3i[k][j][i] = U1d.Bz;
        if (i == ie) pGrid->B1i[k][j][i+1] = Bx;
        if (j == je) pGrid->B2i[k][j+1][i] = U1d.By;
        if (k == ke && ke > ks) pGrid->B3i[k+1][j][i] = U1d.Bz;
        if (x1 > 0.5 && x1 < 1.5) {
          pGrid->B2i[k][j][i] = -U1d.By;
          pGrid->U[k][j][i].B2c = -U1d.By;
          if (j == je) pGrid->B2i[k][j+1][i] = -U1d.By;
        }
      }
    }
  }

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
