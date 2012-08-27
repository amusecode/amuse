#include "copyright.h"
/*============================================================================*/
/*! \file field_loop.c
 *  \brief Problem generator for advection of a field loop test. 
 *
 * PURPOSE: Problem generator for advection of a field loop test.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *   -  problem/rad   = radius of field loop
 *   -  problem/amp   = amplitude of vector potential (and therefore B)
 *   -  problem/vflow = flow velocity
 *   -  problem/drat  = density ratio in loop.  Enables density advection and
 *                      thermal conduction tests.
 *   The flow is automatically set to run along the diagonal. 
 *
 *   Various test cases are possible:
 *   - (iprob=1): field loop in x1-x2 plane (cylinder in 3D)
 *   - (iprob=2): field loop in x2-x3 plane (cylinder in 3D)
 *   - (iprob=3): field loop in x3-x1 plane (cylinder in 3D) 
 *   - (iprob=4): rotated cylindrical field loop in 3D.
 *   - (iprob=5): spherical field loop in rotated plane
 *
 *   A sphere of passive scalar can be added to test advection of scalars.
 *
 *   The temperature in the loop can be changed using drat to test conduction.
 *
 * REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD
 *   via constrined transport", JCP, 205, 509 (2005)			      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3,iprob;
  Real x1c,x2c,x3c,x1f,x2f,x3f;       /* cell- and face-centered coordinates */
  Real x1size,x2size,x3size,lambda=0.0,ang_2=0.0,sin_a2=0.0,cos_a2=1.0,x,y;
  Real rad,amp,vflow,drat,diag;
  Real ***az,***ay,***ax;
#ifdef MHD
  int ku;
#endif
#if (NSCALARS > 0)
  int n;
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  if (((je-js) == 0)) {
    ath_error("[field_loop]: This problem can only be run in 2D or 3D\n");
  }

  if ((ay = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[field_loop]: Error allocating memory for vector pot\n");
  }
  if ((az = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[field_loop]: Error allocating memory for vector pot\n");
  }
  if ((ax = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[field_loop]: Error allocating memory for vector pot\n");
  }

/* Read initial conditions, diffusion coefficients (if needed) */

  rad = par_getd("problem","rad");
  amp = par_getd("problem","amp");
  vflow = par_getd("problem","vflow");
  drat = par_getd_def("problem","drat",1.0);
  iprob = par_getd("problem","iprob");
#ifdef OHMIC
  eta_R = par_getd("problem","eta");
#endif
#ifdef ISOTROPIC_CONDUCTION
  kappa_T = par_getd("problem","kappa");
#endif
#ifdef ANISOTROPIC_CONDUCTION
  chi_C = par_getd("problem","chi");
#endif

/* For (iprob=4) -- rotated cylinder in 3D -- set up rotation angle and
 * wavelength of cylinder */

  if(iprob == 4){
    x1size = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    x3size = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/* We put 1 wavelength in each direction.  Hence the wavelength
 *     lambda = x1size*cos_a;
 *     AND   lambda = x3size*sin_a;
 *     are both satisfied. */

    if(x1size == x3size){
      ang_2 = PI/4.0;
      cos_a2 = sin_a2 = sqrt(0.5);
    }
    else{
      ang_2 = atan(x1size/x3size);
      sin_a2 = sin(ang_2);
      cos_a2 = cos(ang_2);
    }
/* Use the larger angle to determine the wavelength */
    if (cos_a2 >= sin_a2) {
      lambda = x1size*cos_a2;
    } else {
      lambda = x3size*sin_a2;
    }
  }

/* Use vector potential to initialize field loop */

  for (k=ks; k<=ke+1; k++) {
  for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie+1; i++) {
    cc_pos(pGrid,i,j,k,&x1c,&x2c,&x3c);
    x1f = x1c - 0.5*pGrid->dx1;
    x2f = x2c - 0.5*pGrid->dx2;
    x3f = x3c - 0.5*pGrid->dx3;
     
/* (iprob=1): field loop in x1-x2 plane (cylinder in 3D) */

    if(iprob==1) {  
      ax[k][j][i] = 0.0;
      ay[k][j][i] = 0.0;
      if ((x1f*x1f + x2f*x2f) < rad*rad) {
        az[k][j][i] = amp*(rad - sqrt(x1f*x1f + x2f*x2f));
      } else {
        az[k][j][i] = 0.0;
      }
    }

/* (iprob=2): field loop in x2-x3 plane (cylinder in 3D) */

    if(iprob==2) {  
      if ((x2f*x2f + x3f*x3f) < rad*rad) {
        ax[k][j][i] = amp*(rad - sqrt(x2f*x2f + x3f*x3f));
      } else {
        ax[k][j][i] = 0.0;
      }
      ay[k][j][i] = 0.0;
      az[k][j][i] = 0.0;
    }

/* (iprob=3): field loop in x3-x1 plane (cylinder in 3D) */

    if(iprob==3) {  
      if ((x1f*x1f + x3f*x3f) < rad*rad) {
        ay[k][j][i] = amp*(rad - sqrt(x1f*x1f + x3f*x3f));
      } else {
        ay[k][j][i] = 0.0;
      }
      ax[k][j][i] = 0.0;
      az[k][j][i] = 0.0;
    }

/* (iprob=4): rotated cylindrical field loop in 3D.  Similar to iprob=1
 * with a rotation about the x2-axis.  Define coordinate systems (x1,x2,x3)
 * and (x,y,z) with the following transformation rules:
 *    x =  x1*cos(ang_2) + x3*sin(ang_2)
 *    y =  x2
 *    z = -x1*sin(ang_2) + x3*cos(ang_2)
 * This inverts to:
 *    x1  = x*cos(ang_2) - z*sin(ang_2)
 *    x2  = y
 *    x3  = x*sin(ang_2) + z*cos(ang_2)
 */

    if(iprob==4) {
      x = x1c*cos_a2 + x3f*sin_a2;
      y = x2f;
/* shift x back to the domain -0.5*lambda <= x <= 0.5*lambda */
      while(x >  0.5*lambda) x -= lambda;
      while(x < -0.5*lambda) x += lambda;
      if ((x*x + y*y) < rad*rad) {
        ax[k][j][i] = amp*(rad - sqrt(x*x + y*y))*(-sin_a2);
      } else {
        ax[k][j][i] = 0.0;
      }

      ay[k][j][i] = 0.0;

      x = x1f*cos_a2 + x3c*sin_a2;
      y = x2f;
/* shift x back to the domain -0.5*lambda <= x <= 0.5*lambda */
      while(x >  0.5*lambda) x -= lambda;
      while(x < -0.5*lambda) x += lambda;
      if ((x*x + y*y) < rad*rad) {
        az[k][j][i] = amp*(rad - sqrt(x*x + y*y))*(cos_a2);
      } else {
        az[k][j][i] = 0.0;
      }
    }

/* (iprob=5): spherical field loop in rotated plane */

    if(iprob==5) { 
      ax[k][j][i] = 0.0;
      if ((x1f*x1f + x2c*x2c + x3f*x3f) < rad*rad) {
        ay[k][j][i] = amp*(rad - sqrt(x1f*x1f + x2c*x2c + x3f*x3f));
      } else {
        ay[k][j][i] = 0.0;
      }
      if ((x1f*x1f + x2f*x2f + x3c*x3c) < rad*rad) {
        az[k][j][i] = amp*(rad - sqrt(x1f*x1f + x2f*x2f + x3c*x3c));
      } else {
        az[k][j][i] = 0.0;
      }
    }

  }}}

/* Initialize density and momenta.  If drat != 1, then density and temperature
 * will be different inside loop than background values
 */

  x1size = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  x2size = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  x3size = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
  diag = sqrt(x1size*x1size + x2size*x2size + x3size*x3size);
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
     pGrid->U[k][j][i].d = 1.0;
     pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*vflow*x1size/diag;
     pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*vflow*x2size/diag;
     pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*vflow*x3size/diag;

     cc_pos(pGrid,i,j,k,&x1c,&x2c,&x3c);
     if ((x1c*x1c + x2c*x2c + x3c*x3c) < rad*rad) {
       pGrid->U[k][j][i].d = drat;
       pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*vflow*x1size/diag;
       pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*vflow*x2size/diag;
       pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*vflow*x3size/diag;
     }
#if (NSCALARS > 0)
     for (n=0; n<NSCALARS; n++) pGrid->U[k][j][i].s[n] = 0.0;
     if ((x1c*x1c + x2c*x2c + x3c*x3c) < rad*rad) {
       for (n=0; n<NSCALARS; n++)  pGrid->U[k][j][i].s[n] = 1.0;
     }
#endif
  }}}

/* boundary conditions on interface B */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie+1; i++) {
    pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i])/pGrid->dx2 -
                          (ay[k+1][j][i] - ay[k][j][i])/pGrid->dx3;
  }}}
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i])/pGrid->dx3 -
                          (az[k][j][i+1] - az[k][j][i])/pGrid->dx1;
  }}}
  if (ke > ks) {
    ku = ke+1;
  } else {
    ku = ke;
  }
  for (k=ks; k<=ku; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i])/pGrid->dx1 -
                          (ax[k][j+1][i] - ax[k][j][i])/pGrid->dx2;
  }}}
#endif

/* initialize total energy and cell-centered B */

#if defined MHD || !defined ISOTHERMAL
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef MHD
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i  ] + 
                                     pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j  ][i] +
                                     pGrid->B2i[k][j+1][i]);
	if (ke > ks)
	  pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k  ][j][i] + 
                                       pGrid->B3i[k+1][j][i]);
	else
	  pGrid->U[k][j][i].B3c = pGrid->B3i[k][j][i];
#endif

#ifndef ISOTHERMAL
	pGrid->U[k][j][i].E = 1.0/Gamma_1
#ifdef MHD
	  + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
               + SQR(pGrid->U[k][j][i].B3c))
#endif
	  + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
               + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
      }
    }
  }
#endif

  free_3d_array((void***)az);
  free_3d_array((void***)ay);
  free_3d_array((void***)ax);
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
 * current() - computes x3-component of current
 * Bp2()     - computes magnetic pressure (Bx2 + By2)
 * color()   - returns first passively advected scalar s[0]
 * Temperature() - returns temperature for conduction tests
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

#ifdef MHD
/*! \fn static Real current(const GridS *pG, const int i, const int j, const 
 *			   int k)
 *  \brief computes x3-component of current
 */
static Real current(const GridS *pG, const int i, const int j, const int k)
{
  return ((pG->B2i[k][j][i]-pG->B2i[k][j][i-1])/pG->dx1 - 
	  (pG->B1i[k][j][i]-pG->B1i[k][j-1][i])/pG->dx2);
}

/*! \fn static Real Bp2(const GridS *pG, const int i, const int j, const int k)
 *  \brief computes magnetic pressure (Bx2 + By2) */
static Real Bp2(const GridS *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].B1c*pG->U[k][j][i].B1c + 
	  pG->U[k][j][i].B2c*pG->U[k][j][i].B2c);
}

/*! \fn static Real divB(const GridS *pG, const int i, const int j, const int k)
 *  \brief  calculates div(B) */
static Real divB(const GridS *pG, const int i, const int j, const int k)
{
  Real qa;
  if (pG->Nx[2] > 1) {
    qa = (pG->B1i[k][j][i+1]-pG->B1i[k][j][i])/pG->dx1 + 
         (pG->B2i[k][j+1][i]-pG->B2i[k][j][i])/pG->dx2 +
         (pG->B3i[k+1][j][i]-pG->B3i[k][j][i])/pG->dx3;
  } else {
    qa = (pG->B1i[k][j][i+1]-pG->B1i[k][j][i])/pG->dx1 + 
         (pG->B2i[k][j+1][i]-pG->B2i[k][j][i])/pG->dx2;
  }
  return qa;
}
#endif

#if (NSCALARS > 0)
/*! \fn static Real color(const GridS *pG, const int i, const int j,const int k)
 *  \brief returns first passively advected scalar s[0]
static Real color(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0]/pG->U[k][j][i].d;
}
#endif

#ifndef BAROTROPIC
/*! \fn static Real Temperature(const GridS *pG, const int i, const int j, 
 *				const int k)
 *  \brief returns temperature for conduction tests */
static Real Temperature(const GridS *pG, const int i, const int j, const int k)
{
  Real Temp;
  Temp = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*(
    SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3));
#ifdef MHD
  Temp -= 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) 
            + SQR(pG->U[k][j][i].B3c));
#endif
  Temp *= Gamma_1/pG->U[k][j][i].d;
  return Temp;
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#ifdef MHD
  if(strcmp(expr,"J3")==0) return current;
  else if(strcmp(expr,"Bp2")==0) return Bp2;
  else if(strcmp(expr,"divB")==0) return divB;
#endif
#if (NSCALARS > 0)
  if(strcmp(expr,"color")==0) return color;
#endif
  if(strcmp(expr,"Temperature")==0) return Temperature;
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
