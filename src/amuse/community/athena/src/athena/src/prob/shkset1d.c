#include "copyright.h"
/*============================================================================*/
/*! \file shkset1d.c 
 *  \brief Problem generator for 1-D Riemann problems.  
 *
 * PURPOSE: Problem generator for 1-D Riemann problems.  Initial discontinuity
 *   is located so there are equal numbers of cells to the left and right (at
 *   center of grid based on integer index).  Initializes plane-parallel shock
 *   along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:    */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i,il,iu,j,jl,ju,k,kl,ku;
  int shk_dir; /* Shock direction: {1,2,3} -> {x1,x2,x3} */
  Real x1,x2,x3;
  Prim1DS Wl, Wr;
  Cons1DS U1d, Ul, Ur;
  Real Bxl=0.0, Bxr=0.0;

/* Parse left state read from input file: dl,pl,ul,vl,wl,bxl,byl,bzl */

  Wl.d = par_getd("problem","dl");
#ifdef ADIABATIC
  Wl.P = par_getd("problem","pl");
#endif
  Wl.Vx = par_getd("problem","v1l");
  Wl.Vy = par_getd("problem","v2l");
  Wl.Vz = par_getd("problem","v3l");
#ifdef MHD
  Bxl = par_getd("problem","b1l");
  Wl.By = par_getd("problem","b2l");
  Wl.Bz = par_getd("problem","b3l");
#endif
#if (NSCALARS > 0)
  Wl.r[0] = par_getd("problem","r[0]l");
#endif

/* Parse right state read from input file: dr,pr,ur,vr,wr,bxr,byr,bzr */

  Wr.d = par_getd("problem","dr");
#ifdef ADIABATIC
  Wr.P = par_getd("problem","pr");
#endif
  Wr.Vx = par_getd("problem","v1r");
  Wr.Vy = par_getd("problem","v2r");
  Wr.Vz = par_getd("problem","v3r");
#ifdef MHD
  Bxr = par_getd("problem","b1r");
  Wr.By = par_getd("problem","b2r");
  Wr.Bz = par_getd("problem","b3r");
  if (Bxr != Bxl) ath_error(0,"[shkset1d] L/R values of Bx not the same\n");
#endif
#if (NSCALARS > 0)
  Wr.r[0] = par_getd("problem","r[0]r");
#endif

  Ul = Prim1D_to_Cons1D(&Wl, &Bxl);
  Ur = Prim1D_to_Cons1D(&Wr, &Bxr);

/* Parse shock direction */
  shk_dir = par_geti("problem","shk_dir");
  if (shk_dir != 1 && shk_dir != 2 && shk_dir != 3) {
    ath_error("[problem]: shk_dir = %d must be either 1,2 or 3\n",shk_dir);
  }

/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }
  else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
  else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Initialize the grid including the ghost cells.  Discontinuity is always
 * located at x=0, so xmin/xmax in input file must be set appropriately. */

  switch(shk_dir) {
/*--- shock in 1-direction ---------------------------------------------------*/
  case 1:  /* shock in 1-direction  */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive and conserved variables to be L or R state */
          if (x1 <= 0.0) {
            U1d = Ul;
          } else {
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.Mx;
          pGrid->U[k][j][i].M2 = U1d.My;
          pGrid->U[k][j][i].M3 = U1d.Mz;
#ifdef MHD
          pGrid->B1i[k][j][i] = Bxl;
          pGrid->B2i[k][j][i] = U1d.By;
          pGrid->B3i[k][j][i] = U1d.Bz;
          pGrid->U[k][j][i].B1c = Bxl;
          pGrid->U[k][j][i].B2c = U1d.By;
          pGrid->U[k][j][i].B3c = U1d.Bz;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif
        }
      }
    }
    break;

/*--- shock in 2-direction ---------------------------------------------------*/
  case 2:  /* shock in 2-direction  */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive variables to be L or R state */
          if (x2 <= 0.0) {
            U1d = Ul;
          } else {
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.Mz;
          pGrid->U[k][j][i].M2 = U1d.Mx;
          pGrid->U[k][j][i].M3 = U1d.My;
#ifdef MHD
          pGrid->B1i[k][j][i] = U1d.Bz;
          pGrid->B2i[k][j][i] = Bxl;
          pGrid->B3i[k][j][i] = U1d.By;
          pGrid->U[k][j][i].B1c = U1d.Bz;
          pGrid->U[k][j][i].B2c = Bxl;
          pGrid->U[k][j][i].B3c = U1d.By;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif
        }
      }
    }
    break;

/*--- shock in 3-direction ---------------------------------------------------*/
  case 3:  /* shock in 3-direction  */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive variables to be L or R state */
          if (x3 <= 0.0) {
            U1d = Ul;
          } else {
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.My;
          pGrid->U[k][j][i].M2 = U1d.Mz;
          pGrid->U[k][j][i].M3 = U1d.Mx;
#ifdef MHD
          pGrid->B1i[k][j][i] = U1d.By;
          pGrid->B2i[k][j][i] = U1d.Bz;
          pGrid->B3i[k][j][i] = Bxl;
          pGrid->U[k][j][i].B1c = U1d.By;
          pGrid->U[k][j][i].B2c = U1d.Bz;
          pGrid->U[k][j][i].B3c = Bxl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif
        }
      }
    }
  break;
  default:
    ath_error("[shkset1d]: invalid shk_dir = %i\n",shk_dir);
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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}
