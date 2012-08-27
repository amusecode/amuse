#include "copyright.h"
/*============================================================================*/
/*! \file cylfieldloop.c
 *  \brief Problem generator for advection of a field loop test in cylindrical
 *   coordinates. 
 *
 * PURPOSE: Problem generator for advection of a field loop test in cylindrical
 *   coordinates.  Can only be run in 2D or 3D.  Input parameters are:
 *   -  problem/r0     = radial coordinate of loop center
 *   -  problem/phi0   = angular coordinate of loop center
 *   -  problem/rad    = radius of field loop
 *   -  problem/amp    = amplitude of vector potential (and therefore B)
 *   -  problem/omega0 = flow angular velocity
 *   -  problem/vz0    = flow vertical velocity
 *
 * REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD
 *   via constrined transport", JCP, 205, 509 (2005) */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * grav_pot() - gravitational potential
 * grav_acc() - gravitational acceleration
 * M2()       - phi-momentum
 * A3()       - magnetic vector potential
 *============================================================================*/

static Real r0,phi0,amp,rad,rho0,omega0,vz0;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);
Real M2(const Real x1, const Real x2, const Real x3);
Real A3(const Real x1, const Real x2, const Real x3);
static ConsS ***RootSoln=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real x1,x2,x3,lsf,rsf;

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  il = is-nghost;  iu = ie+nghost;
  jl = js-nghost;  ju = je+nghost;
  kl = ks-nghost*(ke>ks);  ku = ke+nghost*(ke>ks);

  nx1 = iu-il+1;
  nx2 = ju-jl+1;
  nx3 = ku-kl+1;

#ifndef MHD
  ath_error("[cylfieldloop]: This problem can only be run in MHD!\n");
#endif

  if ((ie-is)==0 || (je-js)==0) {
    ath_error("[cylfieldloop]: This problem can only be run in 2D or 3D\n");
  }

  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylfieldloop]: Error allocating memory for solution\n");

  /* READ INITIAL CONDITIONS */
  r0     = par_getd("problem","r0");
  phi0   = par_getd("problem","phi0");
  amp    = par_getd("problem","amp");
  rad    = par_getd("problem","rad");
  rho0   = par_getd("problem","rho0");
  omega0 = par_getd("problem","omega0");
  vz0    = par_getd("problem","vz0");

  /* INITIALIZE DENSITY, MOMENTA, INTERFACE FIELDS */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));

        pG->U[k][j][i].d  = rho0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = avg1d(M2,pG,i,j,k);
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz0;
        pG->B1i[k][j][i]  = vecpot2b1i(NULL,A3,pG,i,j,k);
        pG->B2i[k][j][i]  = vecpot2b2i(NULL,A3,pG,i,j,k);
        pG->B3i[k][j][i]  = 0.0;
      }
    }
  }

  /* INITIALIZE TOTAL ENERGY AND CELL-CENTERED FIELDS */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        lsf = (x1 - 0.5*pG->dx1)/x1;
        rsf = (x1 + 0.5*pG->dx1)/x1;

        if (i < iu) 
          pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*pG->B1i[k][j][i+1]);
        else
          pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*vecpot2b1i(NULL,A3,pG,i+1,j,k));

        if (j < ju) 
          pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
        else
          pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + vecpot2b2i(NULL,A3,pG,i,j+1,k));

        if (ke > ks)
          if (k < ku)
            pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
          else
            pG->U[k][j][i].B3c = 0.0;
        else
          pG->U[k][j][i].B3c = pG->B3i[k][j][i];

#ifndef ISOTHERMAL
        pG->U[k][j][i].E = 1.0/Gamma_1 
          + 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c))
          + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#endif /* ISOTHERMAL */

        // SAVE SOLUTION
        RootSoln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  printf("Initial Max divB = %1.10e\n", compute_div_b(pG));

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain, left_x1,  do_nothing_bc);
  bvals_mhd_fun(pDomain, right_x1, do_nothing_bc);

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
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
  printf("Max divB = %1.10e\n", compute_div_b(pM->Domain[0][0].Grid));
}

void Userwork_after_loop(MeshS *pM)
{
  compute_l1_error("CylFieldLoop", pM, RootSoln, 0);
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitatioinal potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega0);
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational acceleration  */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return x1*SQR(omega0);
}

/*! \fn Real M2(const Real x1, const Real x2, const Real x3) 
 *  \brief 2-component of momentum */
Real M2(const Real x1, const Real x2, const Real x3) {
  return rho0*omega0*x1;
}

/*! \fn Real A3(const Real x1, const Real x2, const Real x3) 
 *  \brief 3-component of vector potential */
Real A3(const Real x1, const Real x2, const Real x3) {
  Real X0,X,Y0,Y,dist;

  /* CONVERT TO CARTESIAN COORDINATES */
  X = x1*cos(x2);  Y = x1*sin(x2);

  /* (X0,Y0) IS THE CENTER OF THE LOOP */
  X0 = r0*cos(phi0);  Y0 = r0*sin(phi0);

  /* dist IS DISTANCE TO CENTER OF LOOP */
  dist = sqrt(SQR(X-X0) + SQR(Y-Y0));
  if (dist < rad) {
    return amp*(rad - dist);
  }

  return 0.0;
}
