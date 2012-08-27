#include "copyright.h"
/*============================================================================*/
/*! \file shkset2d.c
 *  \brief Sets up shock at angle to grid to test multidimensional algorithm.
 *
 * PURPOSE: Sets up shock at angle to grid to test multidimensional algorithm.
 *   The grid angle atan(Ly/Lx) is fixed to be atan(0.5), or atan(1), and 
 *   Nx1/Nx2 must be the same ratio as Lx/Ly.  Uses the angle of the shock to
 *   remap ghost cells to the equivalent active grid cells, which requires
 *   that Nx1>32, using special function pointers.  The shock is initialized
 *   with reference to a coordinate system (x,y,z) with transformation rules to
 *   the code coordinate system (x1,x2,x3)
 *   -  x =  x1*cos(alpha) + x2*sin(alpha)
 *   -  y = -x1*sin(alpha) + x2*cos(alpha)
 *   -  z = x3

 *   This inverts to:
 *   -  x1 = x*cos(alpha) - y*sin(alpha)
 *   -  x2 = x*sin(alpha) + y*cos(alpha)
 *   -  x3 = z								  
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - shkset2d_iib() - sets BCs on L-x1 (left edge) of grid.
 * - shkset2d_oib() - sets BCs on R-x1 (right edge) of grid.
 * - shkset2d_ijb() - sets BCs on L-x2 (bottom edge) of grid.
 * - shkset2d_ojb() - sets BCs on R-x2 (top edge) of grid.		      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * shkset2d_iib() - sets BCs on L-x1 (left edge) of grid.
 * shkset2d_oib() - sets BCs on R-x1 (right edge) of grid.
 * shkset2d_ijb() - sets BCs on L-x2 (bottom edge) of grid.
 * shkset2d_ojb() - sets BCs on R-x2 (top edge) of grid.
 *============================================================================*/

void shkset2d_iib(GridS *pGrid);
void shkset2d_oib(GridS *pGrid);
void shkset2d_ijb(GridS *pGrid);
void shkset2d_ojb(GridS *pGrid);

/* Make size of box and dimension of unit cell (r1 x r2) static globals so they
 * can be accessed by boundary value functions */
static Real Lx,Ly;
static int r1,r2;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int kl,ku,irefine,ir,ix1,ix2,nx1,nx2,gcd;
  Real angle, sin_a, cos_a; /* Angle the shock makes with the x1-direction */
  Real rootdx1, rootdx2;
  Prim1DS Wl, Wr;
  Cons1DS Ul, Ur;
  ConsS ql, qr;
  Real Bxl=0.0,Bxr=0.0;
  div_t id;   /* structure containing remainder and quotient */

/* Following are used to compute volume of cell crossed by initial interface
 * that is assigned to left/right states */
  int dll, dlr, drr, drl;
  Real afl_lx, afr_lx, afl_rx, afr_rx;
  Real afl_ly, afr_ly, afl_ry, afr_ry;
  Real vfl, vfr, B1r, B2r;

  if (pGrid->Nx[1] == 1)
    ath_error("[shkset2d]: This problem can only be run in 2D or 3D\n");

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Find number of cells on root grid */

  irefine = 1;
  for (ir=1;ir<=pDomain->Level;ir++) irefine *= 2;

  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  rootdx1 = pGrid->dx1*((double)(irefine)); 
  rootdx2 = pGrid->dx2*((double)(irefine)); 

  nx1 = (int)(Lx/rootdx1);
  nx2 = (int)(Ly/rootdx2);

/* Compute greatest common divisor of nx1,nx2.  The size of the "unit cell"
 * is nx1/gcd by nx2/gcd */

  if((gcd = ath_gcd(nx1,nx2)) < 10)
    ath_error("[shkset2d]: Greatest Common Divisor (nx1,nx2) = %d\n",gcd);

  id = div(nx1,gcd);
  r1 = id.quot;
  if(id.rem != 0)
    ath_error("[shkset2d]: GCD failed, Remainder of %d / %d is %d\n",
	      nx1,gcd,id.rem);

  id = div(nx2,gcd);
  r2 = id.quot;
  if(id.rem != 0)
    ath_error("[shkset2d]: GCD failed, Remainder of %d / %d is %d\n",
	      nx2,gcd,id.rem);

  ath_pout(1,"The unit cell is (%d,1,%d) grid cells in size\n",r1,r2);

/* Compute angle initial interface makes to the grid */

  if(Lx == Ly){
    cos_a = sin_a = sqrt(0.5);
  }
  else{
    angle = atan((double)(Lx/Ly));
    sin_a = sin(angle);
    cos_a = cos(angle);
  }

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
  if (Bxr != Bxl) ath_error(0,"[shkset2d] L/R values of Bx not the same\n");
#endif

  Ul = Prim1D_to_Cons1D(&Wl, &Bxl);
  Ur = Prim1D_to_Cons1D(&Wr, &Bxr);

/* Initialize ql rotated to the (x1,x2,x3) coordinate system */
  ql.d   = Ul.d;
  ql.M1  = Ul.Mx*cos_a - Ul.My*sin_a;
  ql.M2  = Ul.Mx*sin_a + Ul.My*cos_a;
  ql.M3  = Ul.Mz;
#ifdef MHD
  ql.B1c = Bxl*cos_a - Ul.By*sin_a;
  ql.B2c = Bxl*sin_a + Ul.By*cos_a;
  ql.B3c = Ul.Bz;
#endif
#ifdef ADIABATIC
  ql.E   = Ul.E;
#endif

/* Initialize qr rotated to the (x1,x2,x3) coordinate system */
  qr.d   = Ur.d;
  qr.M1  = Ur.Mx*cos_a - Ur.My*sin_a;
  qr.M2  = Ur.Mx*sin_a + Ur.My*cos_a;
  qr.M3  = Ur.Mz;
#ifdef MHD
  qr.B1c = Bxr*cos_a - Ur.By*sin_a;
  qr.B2c = Bxr*sin_a + Ur.By*cos_a;
  qr.B3c = Ur.Bz;
#endif
#ifdef ADIABATIC
  qr.E   = Ur.E;
#endif

/* Initialize the grid */

  for (k=kl; k<=ku; k++) {
    for (j=0; j<=je+nghost; j++) {
      ix2 = j + pGrid->Disp[1];
      for (i=0; i<=ie+nghost; i++) {
	ix1 = i + pGrid->Disp[0];

/* cell is completely in the left state */
	if((drr = r2*(ix1) + r1*(ix2) - gcd*r1*r2) <= 0){
	  pGrid->U[k][j][i] = ql;
#ifdef MHD
	  pGrid->B1i[k][j][i] = ql.B1c;
	  pGrid->B2i[k][j][i] = ql.B2c;
	  pGrid->B3i[k][j][i] = ql.B3c;
#endif /* MHD */
	}
/* cell is completely in the right state */
	else if((dll = r2*(ix1-1) + r1*(ix2-1) - gcd*r1*r2) >= 0){
	  pGrid->U[k][j][i] = qr;
#ifdef MHD
	  pGrid->B1i[k][j][i] = qr.B1c;
	  pGrid->B2i[k][j][i] = qr.B2c;
	  pGrid->B3i[k][j][i] = qr.B3c;
#endif /* MHD */
	}
/* The more complicated case of a cell  split by the interface boundary */
	else{
	  dlr = r2*(ix1-1) + r1*(ix2) - gcd*r1*r2;

	  if(dlr < 0){ /* The boundary hits the right y-face */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = (Real)(-dlr)/(Real)(r2);
	    afr_ry = 1.0 - afl_ry;
	  }
	  else if(dlr > 0){ /* The boundary hits the left x-face */
	    afl_lx = (Real)(-dll)/(Real)(r1);
	    afr_lx = 1.0 - afl_lx;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }
	  else{ /* dlr == 0.0, The boundary hits the grid cell corner */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }

	  drl = r2*(ix1) + r1*(ix2-1) - gcd*r1*r2;

	  if(drl < 0){ /* The boundary hits the right x-face */
	    afl_rx = (Real)(-drl)/(Real)(r1);
	    afr_rx = 1.0 - afl_rx;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }
	  else if(drl > 0){ /* The boundary hits the left y-face */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = (Real)(-dll)/(Real)(r2);
	    afr_ly = 1.0 - afl_ly;
	  }
	  else{ /* drl == 0.0, The boundary hits the grid cell corner */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }

/* The boundary hits both x-interfaces */
	  if(dlr > 0 && drl < 0){ 
	    vfl = 0.5*(afl_lx + afl_rx);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both y-interfaces */
	  else if(dlr < 0 && drl > 0){ 
	    vfl = 0.5*(afl_ly + afl_ry);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both grid cell corners */
	  else if(dlr == 0 && drl == 0){ 
	    vfl = vfr = 0.5;
	  }
/* The boundary hits the left x- and left y-interface */
	  else if(dlr > 0 && drl > 0){
	    vfl = 0.5*afl_lx*afl_ly;
	    vfr = 1.0 - vfl;
	  }
/* dlr<0 && drl<0:  The boundary hits the right x- and right y-interface */
	  else{ 
	    vfr = 0.5*afr_rx*afr_ry;
	    vfl = 1.0 - vfr;
	  }

/* Initialize the x- and y-interface magnetic fields */
#ifdef MHD
	  pGrid->B1i[k][j][i] = afl_lx*ql.B1c + afr_lx*qr.B1c;
	  B1r              = afl_rx*ql.B1c + afr_rx*qr.B1c;

	  pGrid->B2i[k][j][i] = afl_ly*ql.B2c + afr_ly*qr.B2c;
	  B2r              = afl_ry*ql.B2c + afr_ry*qr.B2c;

	  pGrid->B3i[k][j][i] = vfl*ql.B3c + vfr*qr.B3c;
#endif /* MHD */

/* Initialize the volume averaged quantities */
	  pGrid->U[k][j][i].d  = vfl*ql.d + vfr*qr.d;
	  pGrid->U[k][j][i].M1 = vfl*ql.M1 + vfr*qr.M1;
	  pGrid->U[k][j][i].M2 = vfl*ql.M2 + vfr*qr.M2;
	  pGrid->U[k][j][i].M3 = vfl*ql.M3 + vfr*qr.M3;
#ifdef MHD
	  pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i] + B1r);
	  pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i] + B2r);
	  pGrid->U[k][j][i].B3c = vfl*ql.B3c + vfr*qr.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E  = vfl*ql.E + vfr*qr.E;
#endif
	}
      }
    }
  }

/* Set boundary value function pointers */

  bvals_mhd_fun(pDomain, left_x1,  shkset2d_iib);
  bvals_mhd_fun(pDomain, left_x2,  shkset2d_ijb);
  bvals_mhd_fun(pDomain, right_x1, shkset2d_oib);
  bvals_mhd_fun(pDomain, right_x2, shkset2d_ojb);

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

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_iib(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_iib(GridS *pGrid)
{
  const int is = pGrid->is;
  int i, j, k, ju, jl, kl, ku; /* j-upper, j-lower */

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost + r2;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  for (k=kl; k<=ku; k++) {
    for (i=1; i<=nghost; i++) { /* Do NOT Change this loop ordering! */
      for (j=jl; j<=ju; j++) {
	pGrid->U  [k][j][is-i] = pGrid->U  [k][j-r2][is-i+r1];
#ifdef MHD
	pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j-r2][is-i+r1];
	pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j-r2][is-i+r1];
	pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j-r2][is-i+r1];
#endif
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_oib(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_oib(GridS *pGrid)
{
  const int ie = pGrid->ie;
  int i, j, k, ju, jl, kl, ku; /* j-upper, j-lower */

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost - r2;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Note that i=ie+1 is not a boundary condition for the interface field B1i */

  for (k=kl; k<=ku; k++) {
    for (i=1; i<=nghost; i++) { /* Do NOT Change this loop ordering! */
      for (j=jl; j<=ju; j++) {
	pGrid->U[k][j][ie+i] = pGrid->U[k][j+r2][ie+i-r1];
#ifdef MHD
	if(i>1) pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j+r2][ie+i-r1];
	pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j+r2][ie+i-r1];
	pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j+r2][ie+i-r1];
#endif
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_ijb(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_ijb(GridS *pGrid)
{
  const int js = pGrid->js;
  int i, j, k, iu, il, kl, ku; /* i-upper, i-lower */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost + r1;

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  for (k=kl; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
	pGrid->U  [k][js-j][i] = pGrid->U  [k][js-j+r2][i-r1];
#ifdef MHD
	pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js-j+r2][i-r1];
	pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js-j+r2][i-r1];
	pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js-j+r2][i-r1];
#endif
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_ojb(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_ojb(GridS *pGrid)
{
  const int je = pGrid->je;
  int i, j, k, iu, il, kl, ku; /* i-upper, i-lower */

  iu = pGrid->ie + nghost - r1;
  il = pGrid->is - nghost;

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */

  for (k=kl; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
	pGrid->U[k][je+j][i] = pGrid->U[k][je+j-r2][i+r1];
#ifdef MHD
	pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-r2][i+r1];
	if(j>1) pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je+j-r2][i+r1];
	pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je+j-r2][i+r1];
#endif
      }
    }
  }
  return;
}
