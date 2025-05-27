#include "copyright.h"
/*============================================================================*/
/*! \file cylblast.c
 *  \brief Problem generator for blast wave in cylindrical coords.
 *
 * PURPOSE: Problem generator for blast wave in cylindrical coords.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *   -  problem/radius = radius of field initial overpressured region
 *   -  problem/pamb   = ambient pressure
 *   -  problem/prat   = ratio of interior to ambient pressure
 *   -  problem/b0     = initial azimuthal magnetic field (units sqrt(Pamb))
 *   -  problem/rho0   = background density
 *   -  problem/omega0 = initial azimuthal flow angular velocity
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein. */
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
 *============================================================================*/

static Real omega0,rho0;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);
Real M2(const Real x1, const Real x2, const Real x3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real r0,phi0,x0,y0,z0,angle,radius,prat,pamb,b0;
  Real x1,x2,x3,x2i;
  Real x,y,z,Eint,Emag,Ekin;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylblast]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylblast]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylblast]: Only (R,phi) can be used in 2D!\n");
  }

  /* READ IN INITIAL CONDITIONS */
  radius = par_getd("problem","radius");
  pamb   = par_getd("problem","pamb");
  prat   = par_getd("problem","prat");
  rho0   = par_getd("problem","rho0");
  omega0 = par_getd("problem","omega0");
  b0     = par_getd("problem","b0");

  /* PLACEMENT OF CENTER OF BLAST */
  r0   = par_getd("problem","r0");
  phi0 = par_getd("problem","phi0");
  z0   = par_getd("problem","z0");

  /* ORIENTATION OF FIELD W.R.T. POS. X-AXIS */
  angle = (PI/180.0)*par_getd("problem","angle");

  x0 = r0*cos(phi0);
  y0 = r0*sin(phi0);

  /* SET UP UNIFORM AMBIENT MEDIUM WITH CIRCULAR OVER-PRESSURED REGION */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        x2i = x2 - 0.5*pG->dx2;

        pG->U[k][j][i].d  = rho0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega0;
        pG->U[k][j][i].M3 = 0.0;
#ifdef MHD
        /* SET UP A PLANAR MAGNETIC FIELD IN THE X-Y (R-PHI) PLANE */
        pG->B1i[k][j][i]   = b0*(cos(angle)*cos(x2)+sin(angle)*sin(x2)); 
        pG->B2i[k][j][i]   = b0*(-cos(angle)*sin(x2i)+sin(angle)*cos(x2i));
        pG->B3i[k][j][i]   = 0.0;
        pG->U[k][j][i].B1c = b0*(cos(angle)*cos(x2)+sin(angle)*sin(x2));
        pG->U[k][j][i].B2c = b0*(-cos(angle)*sin(x2)+sin(angle)*cos(x2));
        pG->U[k][j][i].B3c = 0.0;
#endif
        /* CARTESIAN POSITION OF CELL CENTER */
        x = x1*cos(x2);
        y = x1*sin(x2);
        z = x3;

        /* INITIALIZE TOTAL ENERGY */
#ifndef ISOTHERMAL
        Eint = pamb/Gamma_1;
        if (SQR(x-x0) + SQR(y-y0) + SQR(z-z0) < SQR(radius)) {
          Eint *= prat;
        }
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
        pG->U[k][j][i].E = Eint + Emag + Ekin;
#endif /* ISOTHERMAL */
      }
    }
  }

  /* Enroll the gravitational function and radial BC */
  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain,left_x1, do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,left_x2, do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x2,do_nothing_bc);

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
}

void Userwork_after_loop(MeshS *pM)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) {
 *  \brief  Gravitational potential*/
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega0);
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3) {
 *  \brief Gravitational acceleration  */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return x1*SQR(omega0);
}

/*! \fn Real M2(const Real x1, const Real x2, const Real x3) 
 *  \brief 2-component of momentum */
Real M2(const Real x1, const Real x2, const Real x3) {
  return rho0*omega0*x1;
}
