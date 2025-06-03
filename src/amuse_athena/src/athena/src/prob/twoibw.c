#include "copyright.h"
/*============================================================================*/
/*! \file twoibw.c
 *  \brief Problem generator for two interacting blast waves test.
 *
 * PURPOSE: Problem generator for two interacting blast waves test.
 *
 * REFERENCE: P. Woodward & P. Colella, "The numerical simulation of 
 *   two-dimensional fluid flow with strong shocks", JCP, 54, 115, sect. IVa  */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef ISOTHERMAL
#error : The twoibw test only works for adiabatic EOS.
#endif /* ISOTHERMAL */
#ifndef HYDRO
#error : The twoibw test only works for hydro.
#endif /* HYDRO */

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int shk_dir;
  Real x1,x2,x3;

  shk_dir = par_geti("problem","shk_dir");
  if (shk_dir < 1 || shk_dir > 3) {
    ath_error("[problem]: shk_dir = %d must be either 1, 2 or 3\n",shk_dir);
  }

  if(shk_dir == 1) {
/* setup dependent variables in X1-direction */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
	for (i=is; i<=ie; i++) {
/* Calculate the cell center position of the cell i,j,k */
	  cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	  pGrid->U[k][j][i].d = 1.0;
	  pGrid->U[k][j][i].M1 = 0.0;
	  pGrid->U[k][j][i].M2 = 0.0;
	  pGrid->U[k][j][i].M3 = 0.0;
	  if (x1 < 0.1) {
	    pGrid->U[k][j][i].E = 1.0e3/Gamma_1 ;
	  }
	  else if (x1 > 0.9) {
	    pGrid->U[k][j][i].E = 1.0e2/Gamma_1 ;
	  }
	  else {
	    pGrid->U[k][j][i].E = 0.01/Gamma_1 ;
	  }
	}
      }
    }
  }
  else if (shk_dir == 2) {
/* setup dependent variables in X2-direction */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
	for (i=is; i<=ie; i++) {
/* Calculate the cell center position of the cell i,j,k */
	  cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	  pGrid->U[k][j][i].d = 1.0;
	  pGrid->U[k][j][i].M1 = 0.0;
	  pGrid->U[k][j][i].M2 = 0.0;
	  pGrid->U[k][j][i].M3 = 0.0;
	  if (x2 < 0.1) {
	    pGrid->U[k][j][i].E = 1.0e3/Gamma_1 ;
	  }
	  else if (x2 > 0.9) {
	    pGrid->U[k][j][i].E = 1.0e2/Gamma_1 ;
	  }
	  else {
	    pGrid->U[k][j][i].E = 0.01/Gamma_1 ;
	  }
	}
      }
    }
  }
  else { /* shk_dir == 3 */
/* setup dependent variables in X3-direction */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
	for (i=is; i<=ie; i++) {
/* Calculate the cell center position of the cell i,j,k */
	  cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	  pGrid->U[k][j][i].d = 1.0;
	  pGrid->U[k][j][i].M1 = 0.0;
	  pGrid->U[k][j][i].M2 = 0.0;
	  pGrid->U[k][j][i].M3 = 0.0;
	  if (x3 < 0.1) {
	    pGrid->U[k][j][i].E = 1.0e3/Gamma_1 ;
	  }
	  else if (x3 > 0.9) {
	    pGrid->U[k][j][i].E = 1.0e2/Gamma_1 ;
	  }
	  else {
	    pGrid->U[k][j][i].E = 0.01/Gamma_1 ;
	  }
	}
      }
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
