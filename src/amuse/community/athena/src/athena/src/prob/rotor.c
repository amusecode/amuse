#include "copyright.h"
/*============================================================================*/
/*! \file rotor.c
 *  \brief Sets up 2D rotor test problem.
 *
 * PURPOSE: Sets up 2D rotor test problem.  The center of the grid is assumed to
 *   have coordinates (x1,x2) = [0,0]; the grid initialization must be
 *   consistent with this
 *
 * REFERENCE: G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes",
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
#error : The rotor problem can only be run with MHD.
#endif
#ifdef ISOTHERMAL 
#error : The rotor problem can only be run with an ADIABATIC eos.
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,is,ie,js,je,ks,ke;
  Real v0,p0,bx0,x1,x2,x3,rad,frac,r0,r1;

/* Read initial conditions from 'athinput' */

  v0 = par_getd("problem","v0");
  p0 = par_getd("problem","p0");
  bx0 = par_getd("problem","bx0");
  r0 = par_getd("problem","r0");
  r1 = par_getd("problem","r1");

/* Initialize the grid.  Note the center is always assumed to have coordinates
 * x1=0, x2=0; the grid range in the input file must be consistent with this */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d = 1.0;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
        pGrid->B1i[k][j][i] = bx0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->U[k][j][i].B1c = bx0;
        pGrid->U[k][j][i].B2c = 0.0;
        pGrid->U[k][j][i].B3c = 0.0;

/* reset density, velocity if cell is inside rotor */

        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2);
        if (rad <= r0) {
          pGrid->U[k][j][i].d  = 10.0;
          pGrid->U[k][j][i].M1 = -100.0*v0*x2;
          pGrid->U[k][j][i].M2 = 100.0*v0*x1;
        } else {

/* smooth solution between r0 and r1.  For no smoothing, set r1<0 in input */

          if (rad <= r1) {
            frac = (0.115 - rad)/(0.015);
            pGrid->U[k][j][i].d  = 1.0 + 9.0*frac;
            pGrid->U[k][j][i].M1 = -frac*100.0*v0*x2;
            pGrid->U[k][j][i].M2 = frac*100.0*v0*x1;
          }
        }

        pGrid->U[k][j][i].E = p0/Gamma_1 + 0.5*bx0*bx0
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2))
          /pGrid->U[k][j][i].d;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] = bx0;
    }
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
