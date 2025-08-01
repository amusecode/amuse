#include "copyright.h"
/*============================================================================*/
/*! \file lw_implode.c
 *  \brief Problem generator for square implosion problem.
 *
 * PURPOSE: Problem generator for square implosion problem.
 *
 * REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)    */
/*============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MHD
#error : The lw_implode problem only works for hydro.
#endif

/* computes difference d{i,j}-d{j,i} to test if solution is symmetric */
static Real expr_diff_d(const GridS *pG, const int i, const int j, const int k);

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int ir,irefine,nx2;
  Real d_in,p_in,d_out,p_out,Ly,rootdx2;

/* Set up the grid bounds for initializing the grid */
  if (pGrid->Nx[0] <= 1 || pGrid->Nx[1] <= 1) {
    ath_error("[problem]: This problem requires Nx1 > 1, Nx2 > 1\n");
  }

  d_in = par_getd("problem","d_in");
  p_in = par_getd("problem","p_in");

  d_out = par_getd("problem","d_out");
  p_out = par_getd("problem","p_out");

/* Find number of Nx2 cells on root grid.  At x=0, interface is at nx2/2 */

  irefine = 1;
  for (ir=1;ir<=pDomain->Level;ir++) irefine *= 2;

  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  rootdx2 = pGrid->dx2*((double)(irefine));
  nx2 = (int)(Ly/rootdx2);
  nx2 /= 2;

/* Initialize the grid */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;

	if(((j-js + pDomain->Disp[1])+(i-is + pDomain->Disp[0])) > (nx2*irefine)) {
	  pGrid->U[k][j][i].d  = d_out;
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E = p_out/Gamma_1;
#endif
	} else {
	  pGrid->U[k][j][i].d  = d_in;
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E = p_in/Gamma_1;
#endif
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
  if(strcmp(expr,"diff_d")==0) return expr_diff_d;
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
/*! \fn static Real expr_diff_d(const GridS *pG, const int i, const int j, 
 *				const int k)
 *  \brief computes difference d{i,j}-d{j,i} to test if solution is symmetric */
static Real expr_diff_d(const GridS *pG, const int i, const int j, const int k)
{
 return (pG->U[k][j][i].d - pG->U[k][i][j].d);
}
