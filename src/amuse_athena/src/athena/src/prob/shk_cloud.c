#include "copyright.h"
/*============================================================================*/
/*! \file shk_cloud.c
 *  \brief Problem generator for shock-cloud problem; planar shock impacting
 *   a single spherical cloud.
 *
 * PURPOSE: Problem generator for shock-cloud problem; planar shock impacting
 *   a single spherical cloud.  Input parameters are:
 *    - problem/Mach   = Mach number of incident shock
 *    - problem/drat   = density ratio of cloud to ambient
 *    - problem/beta   = ratio of Pgas/Pmag
 *    - problem/iprob  = integer flag to determine problem
 *
 *   The cloud radius is fixed at 1.0.  The center of the coordinate system
 *   defines the center of the cloud, and should be in the middle of the cloud.
 *   The shock is initially at x1=-2.0.  A typical grid domain should span
 *   x1 in [-3.0,7.0] , y and z in [-2.5,2.5] (see input file in /tst)
 *   Various test cases are possible:
 *   - (iprob=1): B parallel to shock normal
 *   - (iprob=2): B perpendicular to shock normal -- NOT YET IMPLEMENTED
 *
 *   If the code is configured with nscalars>0, the cloud material is labeled
 *   with U[k][j][i].s[0]=1.						      
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * shk_cloud_iib() - fixes BCs on L-x1 (left edge) of grid to postshock flow. */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* postshock flow variables are shared with IIB function */

static Real dl,pl,ul;
#ifdef MHD
static Real bxl,byl,bzl;
#endif /* MHD */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * shk_cloud_iib() - fixes BCs on L-x1 (left edge) of grid to postshock flow.
 *============================================================================*/

void shk_cloud_iib(GridS *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  Real jump1, jump2, jump3;
  Real x1,x2,x3,diag,rad,xshock;
  Real Mach,drat;
#ifdef MHD
  Real beta,bxr,byr,bzr;
#endif /* MHD */
  Real dr,pr,ur;

/* Read input parameters */

  xshock = -2.0;
  rad    = 1.0;
  Mach = par_getd("problem","Mach");
  drat = par_getd("problem","drat");
  iprob = par_geti("problem","iprob");
  iprob = 1;
#ifdef MHD
  beta = par_getd("problem","beta");
#endif
  
/* Set paramters in ambient medium ("R-state" for shock) */

  dr = 1.0;
  pr = 1.0/Gamma;
  ur = 0.0;

/* Uses Rankine Hugoniot relations for adiabatic gas to initialize problem */

  jump1 = (Gamma + 1.0)/(Gamma_1 + 2.0/(Mach*Mach));
  jump2 = (2.0*Gamma*Mach*Mach - Gamma_1)/(Gamma + 1.0);
  jump3 = 2.0*(1.0 - 1.0/(Mach*Mach))/(Gamma + 1.0);

  dl = dr*jump1;
  pl = pr*jump2;
#ifdef ISOTHERMAL
  ul = ur + jump3*Mach*Iso_csound;
#else
  ul = ur + jump3*Mach*sqrt(Gamma*pr/dr);
#endif

/* Initialize the grid */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        diag = sqrt(x1*x1 + x2*x2 + x3*x3);

/* postshock flow */
        if(x1 < xshock) {
          pGrid->U[k][j][i].d  = dl;
          pGrid->U[k][j][i].M1 = ul*dl;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef MHD
          pGrid->B1i[k][j][i] = bxl;
          pGrid->B2i[k][j][i] = byl;
          pGrid->B3i[k][j][i] = bzl;
          pGrid->U[k][j][i].B1c = bxl;
          pGrid->U[k][j][i].B2c = byl;
          pGrid->U[k][j][i].B3c = bzl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pl/Gamma_1 
#ifdef MHD
            + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
            + 0.5*dl*(ul*ul);
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = 0.0;
#endif

/* preshock ambient gas */
        } else {
          pGrid->U[k][j][i].d  = dr;
          pGrid->U[k][j][i].M1 = ur*dr;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef MHD
          pGrid->B1i[k][j][i] = bxr;
          pGrid->B2i[k][j][i] = byr;
          pGrid->B3i[k][j][i] = bzr;
          pGrid->U[k][j][i].B1c = bxr;
          pGrid->U[k][j][i].B2c = byr;
          pGrid->U[k][j][i].B3c = bzr;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pr/Gamma_1
#ifdef MHD
            + 0.5*(bxr*bxr + byr*byr + bzr*bzr)
#endif
            + 0.5*dr*(ur*ur);
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = 0.0;
#endif
        }

/* cloud interior */
        if (diag < rad) {
          pGrid->U[k][j][i].d  = dr*drat;
          pGrid->U[k][j][i].M1 = ur*dr*drat;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef MHD
          pGrid->B1i[k][j][i] = bxr;
          pGrid->B2i[k][j][i] = byr;
          pGrid->B3i[k][j][i] = bzr;
          pGrid->U[k][j][i].B1c = bxr;
          pGrid->U[k][j][i].B2c = byr;
          pGrid->U[k][j][i].B3c = bzr;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pr/Gamma_1
#ifdef MHD
            + 0.5*(bxr*bxr + byr*byr + bzr*bzr)
#endif
            + 0.5*dr*drat*(ur*ur);
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = drat;
#endif
        }
      }
    }
  }

/* boundary conditions on interface B */

#ifdef MHD
  i = ie+1;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][i] = bxr;
    }
  }
  j = je+1;
  for (k=ks; k<=ke; k++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      if (x1 < xshock) {pGrid->B2i[k][j][i] = byl;}
      else {pGrid->B2i[k][j][i] = byr;}
    }
  }
  if (ke > ks) {
    k = ke+1;
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        if (x1 < xshock) {pGrid->B3i[k][j][i] = bzl;}
        else {pGrid->B3i[k][j][i] = bzr;}
      }
    }
  }
#endif

/* Set IIB value function pointer */

  if (pDomain->Disp[0] == 0) bvals_mhd_fun(pDomain,left_x1,shk_cloud_iib);

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
 * color()   - returns first passively advected scalar s[0]
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

#if (NSCALARS > 0)
static Real color(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"color")==0) return color;
#endif

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

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void shk_cloud_iib(GridS *pGrid)
 *  \brief Sets boundary condition on left X boundary (iib) 
 *
 * Note quantities at this boundary are held fixed at the downstream state
 */

void shk_cloud_iib(GridS *pGrid)
{
  int i=0,j=0,k=0;
  int js,je,ks,ke;

  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][i].d  = dl;
        pGrid->U[k][j][i].M1 = ul*dl;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
#ifdef MHD
        pGrid->B1i[k][j][i] = bxl;
        pGrid->B2i[k][j][i] = byl;
        pGrid->B3i[k][j][i] = bzl;
        pGrid->U[k][j][i].B1c = bxl;
        pGrid->U[k][j][i].B2c = byl;
        pGrid->U[k][j][i].B3c = bzl;
#endif
#ifdef ADIABATIC
        pGrid->U[k][j][i].E = pl/Gamma_1
#ifdef MHD
          + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
          + 0.5*dl*(ul*ul);
#endif
#if (NSCALARS > 0)
        pGrid->U[k][j][i].s[0] = 0.0;
#endif
      }
    }
  }
}
