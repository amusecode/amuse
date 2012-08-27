#include "copyright.h"
/*============================================================================*/
/*! \file shu-osher.c
 *  \brief Problem generator for Shu-Osher shocktube test, involving
 *   interaction of a Mach 3 shock with a sine wave density distribution.  
 *
 * PURPOSE: Problem generator for Shu-Osher shocktube test, involving
 *   interaction of a Mach 3 shock with a sine wave density distribution.  
 *
 * REFERENCE: C.W. Shu & S. Osher, "Efficient implementation of essentially
 *   non-oscillatory shock-capturing schemes, II", JCP, 83, 32 (1998)	      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef ISOTHERMAL
#error : The shu-osher test only works for adiabatic EOS.
#endif /* ISOTHERMAL */
#ifndef HYDRO
#error : The shu-osher test only works for hydro.
#endif /* HYDRO */

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i,il,iu,js,ks;
  Real dl,pl,ul,vl,wl,x1,x2,x3;

/* Set up the grid bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  js = pGrid->js;
  ks = pGrid->ks;

  if (pGrid->Nx[1] > 1 || pGrid->Nx[2] > 1) {
    ath_error("Shu Osher test only works for 1D problem in x1-direction\n");
  }

/* setup dependent variables */
  dl = 3.857143;
  ul = 2.629369;
  pl = 10.33333;
  wl = 0.0;
  vl = 0.0;

  for (i=il; i<=iu; i++) {
/* Calculate the cell center position of the cell i,j */
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);

    if (x1 < -0.8) {
      pGrid->U[ks][js][i].d = dl;
      pGrid->U[ks][js][i].M1 = ul*dl;
      pGrid->U[ks][js][i].M2 = vl*dl;
      pGrid->U[ks][js][i].M3 = wl*dl;
      pGrid->U[ks][js][i].E = pl/Gamma_1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
    }
    else {
      pGrid->U[ks][js][i].d = 1.0 + 0.2*sin(5.0*PI*x1);
      pGrid->U[ks][js][i].M1 = 0.0;
      pGrid->U[ks][js][i].M2 = 0.0;
      pGrid->U[ks][js][i].M3 = 0.0;
      pGrid->U[ks][js][i].E = 1.0/Gamma_1;
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
