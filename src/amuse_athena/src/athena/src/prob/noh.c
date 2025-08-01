#include "copyright.h"
/*============================================================================*/
/*! \file noh.c
 *  \brief Spherical Noh implosion problem, from Liska & Wendroff, section 4.5
 *   (figure 4.7).
 *
 * PURPOSE: Spherical Noh implosion problem, from Liska & Wendroff, section 4.5
 *   (figure 4.7).  Tests code on VERY strong shock, also sensitive to 
 *   carbuncle instability.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - void noh3d_oib() - sets BCs on R-x1 boundary
 * - void noh3d_ojb() - sets BCs on R-x2 boundary
 * - void noh3d_okb() - sets BCs on R-x3 boundary
 * - void scat_plot() - makes scatter plot of density
 *
 * REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)    */
/*============================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * void noh3d_oib() - sets BCs on R-x1 boundary
 * void noh3d_ojb() - sets BCs on R-x2 boundary
 * void noh3d_okb() - sets BCs on R-x3 boundary
 * void scat_plot() - makes scatter plot of density
 *============================================================================*/

void noh3d_oib(GridS *pGrid);
void noh3d_ojb(GridS *pGrid);
void noh3d_okb(GridS *pGrid);

#ifdef MHD
#error : This is not a MHD problem.
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pD)
{
  GridS *pGrid=(pD->Grid);
  int i,j,k;
  Real x1,x2,x3,r;

/* Initialize the grid: d=1, v=-1.0 in radial direction, p=10^-6 */

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	if (pGrid->Nx[2] > 1) {
          r = sqrt(x1*x1 + x2*x2 + x3*x3);
        } else {
          r = sqrt(x1*x1 + x2*x2);
        }
	pGrid->U[k][j][i].d = 1.0;
	pGrid->U[k][j][i].M1 = -x1/r;
	pGrid->U[k][j][i].M2 = -x2/r;
	if (pGrid->Nx[2] > 1) {
	  pGrid->U[k][j][i].M3 = -x3/r;
        } else {
	  pGrid->U[k][j][i].M3 = 0.0;
        }
	pGrid->U[k][j][i].E = 1.0e-6/Gamma_1 + 0.5;
      }
    }
  }

/* Set boundary value function pointers */

  bvals_mhd_fun(pD, right_x1,noh3d_oib);
  bvals_mhd_fun(pD, right_x2,noh3d_ojb);
  if (pGrid->Nx[2] > 1) bvals_mhd_fun(pD, right_x3,noh3d_okb);

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

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void noh3d_oib(GridS *pGrid)
 *  \brief Sets boundary condition on right X1 boundary (oib) for noh3d test
 *
 * Note quantities at this boundary are held fixed at the time-dependent
 * upstream state
 */

void noh3d_oib(GridS *pGrid)
{
  int i, ie = pGrid->ie;
  int j, jl = pGrid->js - nghost, ju = pGrid->je + nghost;
  int k, kl, ku;
  Real x1,x2,x3,r,d0,f_t;

  if (pGrid->Nx[2] > 1) {
    kl = pGrid->ks - nghost;
    ku = pGrid->ke + nghost;
  } else {
    kl = pGrid->ks;
    ku = pGrid->ke;
  }

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju;  j++) {
      for (i=ie+1;  i<=ie+nghost;  i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        if (pGrid->Nx[2] > 1) {
	  r = sqrt(x1*x1 + x2*x2 + x3*x3);
          f_t = (1.0 + pGrid->time/r)*(1.0 + pGrid->time/r);
        } else {
	  r = sqrt(x1*x1 + x2*x2);
          f_t = (1.0 + pGrid->time/r);
        }
	d0 = 1.0*f_t;
   
	pGrid->U[k][j][i].d  = d0;
	pGrid->U[k][j][i].M1 = -x1*d0/r;
	pGrid->U[k][j][i].M2 = -x2*d0/r;
        if (pGrid->Nx[2] > 1) {
	  pGrid->U[k][j][i].M3 = -x3*d0/r;
          pGrid->U[k][j][i].E = 1.0e-6*pow(f_t,1.0+Gamma)/Gamma_1 + 0.5*d0;
        } else {
	  pGrid->U[k][j][i].M3 = 0.0;
          pGrid->U[k][j][i].E = 1.0e-6/Gamma_1 + 0.5*d0;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void noh3d_ojb(GridS *pGrid)
 *  \brief Sets boundary condition on right X2 boundary (ojb) for noh3d test
 * 
 * Note quantities at this boundary are held fixed at the time-dependent
 * upstream state
 */

void noh3d_ojb(GridS *pGrid)
{
  int j, je = pGrid->je;
  int i, il = pGrid->is - nghost, iu = pGrid->ie + nghost;
  int k, kl, ku;
  Real x1,x2,x3,r,d0,f_t;

  if (pGrid->Nx[2] > 1) {
    kl = pGrid->ks - nghost;
    ku = pGrid->ke + nghost;
  } else {
    kl = pGrid->ks;
    ku = pGrid->ke;
  }

  for (k=kl; k<=ku; k++) {
    for (j=je+1; j<=je+nghost; j++) {
      for (i=il; i<=iu; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        if (pGrid->Nx[2] > 1) {
	  r = sqrt(x1*x1 + x2*x2 + x3*x3);
          f_t = (1.0 + pGrid->time/r)*(1.0 + pGrid->time/r);
        } else {
	  r = sqrt(x1*x1 + x2*x2);
          f_t = (1.0 + pGrid->time/r);
        }
	d0 = 1.0*f_t;

	pGrid->U[k][j][i].d  = d0;
	pGrid->U[k][j][i].M1 = -x1*d0/r;
	pGrid->U[k][j][i].M2 = -x2*d0/r;
        if (pGrid->Nx[2] > 1) {
	  pGrid->U[k][j][i].M3 = -x3*d0/r;
          pGrid->U[k][j][i].E = 1.0e-6*pow(f_t,1.0+Gamma)/Gamma_1 + 0.5*d0;
        } else {
	  pGrid->U[k][j][i].M3 = 0.0;
          pGrid->U[k][j][i].E = 1.0e-6/Gamma_1 + 0.5*d0;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*! \fn void noh3d_okb(GridS *pGrid) 
 *  \brief Sets boundary condition on right X3 boundary (okb) for noh3d test
 * 
 * Note quantities at this boundary are held fixed at the time-dependent
 * upstream state
 */

void noh3d_okb(GridS *pGrid)
{
  int i, il = pGrid->is - nghost, iu = pGrid->ie + nghost;
  int j, jl = pGrid->js - nghost, ju = pGrid->je + nghost;
  int k, ke = pGrid->ke;
  Real x1,x2,x3,r,d0,f_t;

  for (k=ke+1; k<=ke+nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	r = sqrt(x1*x1 + x2*x2 + x3*x3);
	f_t = (1.0 + pGrid->time/r)*(1.0 + pGrid->time/r);
	d0 = 1.0*f_t;

	pGrid->U[k][j][i].d  = d0;
	pGrid->U[k][j][i].M1 = -x1*d0/r;
	pGrid->U[k][j][i].M2 = -x2*d0/r;
	pGrid->U[k][j][i].M3 = -x3*d0/r;
	pGrid->U[k][j][i].E = 1.0e-6*pow(f_t,1.0+Gamma)/Gamma_1 + 0.5*d0;
      }
    }
  }
}
